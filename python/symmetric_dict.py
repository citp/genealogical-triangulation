class SymmetricDict(dict):

    def __getitem__(self, key):
        return dict.__getitem__(self,
                                key if key[0] < key[1] else (key[1], key[0]))

    def __setitem__(self, key, value):
        dict.__setitem__(self, key if key[0] < key[1] else (key[1],key[0]),
                         value)

    def __delitem__(self, key):
        dict.__delitem__(self, key if key[0] < key[1] else (key[1],key[0]))

    def __contains__(self, key):
        return dict.__contains__(self, key if key[0] < key[1] else (key[1],key[0]))
