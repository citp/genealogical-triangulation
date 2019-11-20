from os import SEEK_END, SEEK_CUR, SEEK_SET
from struct import iter_unpack, unpack
from collections import defaultdict

import numpy as np

FLOAT_SIZE = 8

class BinarySimulationDeserializer:
    def __init__(self, filename):
        self._simulation_file = open(filename, "rb")
        self._end = _last_byte(self._simulation_file)
        header, header_length = _read_header(self._simulation_file)
        self._header_length = header_length
        starts, unlabeled_order = _offsets(header)
        self._starts = starts
        self._unlabeled_order = unlabeled_order
        self._num_elements = sum(len(x)
                                 for x
                                 in self._unlabeled_order.values())
        buf_size = max(len(x) for x
                       in self._unlabeled_order.values()) * FLOAT_SIZE
        self._simulation_file = open(filename, "rb", buffering = buf_size)
        assert len(header) == self._num_elements
        assert ((self._end - header_length) / 8) % self._num_elements == 0

    @property
    def anchors(self):
        return set(self._starts.keys())

    def anchor_shared(self, anchor):
        shared_data = np.array(self.read_all_for_anchor(anchor))
        unlabeled_order = self._unlabeled_order[anchor]
        num_unlabeled = len(unlabeled_order)
        ret = dict()
        data_len = len(shared_data)
        for i, unlabeled in enumerate(unlabeled_order):
            ret[unlabeled] = shared_data[i:data_len:num_unlabeled]
        return ret

    def read_all_for_anchor(self, anchor):
        assert anchor in self._starts.keys()
        start = self._starts[anchor]
        unlabeled_order = self._unlabeled_order[anchor]
        num_unlabeled = len(unlabeled_order)
        read_length = num_unlabeled * FLOAT_SIZE
        skip_length = (self._num_elements - num_unlabeled) * FLOAT_SIZE
        self._simulation_file.seek(self._header_length + start * FLOAT_SIZE)
        data = []
        unpack_str = "d" * num_unlabeled
        while self._simulation_file.tell() < self._end:
            #assert self._simulation_file.tell() + read_length <= self._end
            shared_bytes = self._simulation_file.read(read_length)
            self._simulation_file.seek(skip_length, SEEK_CUR)
            floats = unpack(unpack_str, shared_bytes)
            data.extend(floats)
        return data

    def anchor_shared_alt(self, anchor):
        shared_dict = defaultdict(list)
        for unlabeled, shared in self.iter_shared_for_anchor(anchor):
            shared_dict[unlabeled].append(shared)
        return dict(shared_dict)

    def iter_shared_for_anchor(self, anchor):
        assert anchor in self._starts.keys()
        start = self._starts[anchor]
        unlabeled_order = self._unlabeled_order[anchor]
        num_unlabeled = len(unlabeled_order)
        read_length = num_unlabeled * FLOAT_SIZE
        skip_length = (self._num_elements - num_unlabeled) * FLOAT_SIZE
        self._simulation_file.seek(self._header_length + start * FLOAT_SIZE)
        while self._simulation_file.tell() < self._end:
            #assert self._simulation_file.tell() + read_length <= self._end
            shared_bytes = self._simulation_file.read(read_length)
            self._simulation_file.seek(skip_length, SEEK_CUR)
            floats = iter_unpack("d", shared_bytes)
            for unlabeled, shared in zip(unlabeled_order, floats):
                yield (unlabeled, shared[0])
        
        
def _last_byte(file_object):
    file_object.seek(0, SEEK_END)
    end = file_object.tell()
    file_object.seek(0)
    return end

def _offsets(header):
    anchor_starts = dict()
    anchor_unlabeled = dict()
    current_anchor = None
    current_offset = 0
    current_unlabeled = []
    for anchor, unlabeled in header:
        if current_anchor is None:
            current_anchor = anchor
        if anchor != current_anchor:
            anchor_starts[current_anchor] = current_offset
            anchor_unlabeled[current_anchor] = current_unlabeled
            current_offset += len(current_unlabeled)
            current_unlabeled = []
            current_anchor = anchor
        current_unlabeled.append(unlabeled)
    anchor_starts[current_anchor] = current_offset
    anchor_unlabeled[current_anchor] = current_unlabeled
    assert len(header) == current_offset + len(current_unlabeled)
    return (anchor_starts, anchor_unlabeled)

def _read_header(file_object):
    file_object.seek(0, SEEK_SET)
    # Length of header in bytes. This does not include the size field itself.
    header_length = int.from_bytes(file_object.read(8), byteorder = "little")
    header_bytes = file_object.read(header_length)
    pairs = list(iter_unpack("<LL", header_bytes))
    return (pairs, file_object.tell())
