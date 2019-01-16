from random import choice

from sex import Sex, SEXES

STR_BASE = "<Node id = {}, mother = {}, father = {}"

class NodeGenerator:
    """
    We use a node generator so that nodes point to each other by id
    rather than by reference. Pickle does not handle highly recursive
    datastructures well. This is a performance optimization that
    should allow pickle and other modules that use recursion to work
    on nodes and populations without overflowing the stack.  The IDs
    are also usefull when using multiprocessing, as just the Node's id
    can be passed from one process to the other.
    """
    def __init__(self):
        self._id = 0
        self._mapping = dict()
        
    def generate_node(self, father = None, mother = None,
                      suspected_father = None, suspected_mother = None,
                      sex = None, twin = None):
        node_id = self._id
        self._id += 1
        if father is not None:
            father_id = father._id
        else:
            father_id = None
        if suspected_father is not None:
            suspected_father_id = suspected_father._id
        elif father is not None:
            suspected_father_id = father._id
        else:
            suspected_father_id = None
            
        if mother is not None:
            mother_id = mother._id
        else:
            mother_id = None
        if suspected_mother is not None:
            suspected_mother_id = suspected_mother._id
        elif mother is not None:
            suspected_mother_id = mother._id
        else:
            suspected_mother_id = None

        if twin is not None:
            twin_id = twin._id
        else:
            twin_id = None

        node = Node(self, node_id, father_id, mother_id, suspected_father_id,
                    suspected_mother_id, sex = sex, twin_id = twin_id)
        self._mapping[node_id] = node
        return node

    def twin_node(self, template):
        mother = template.mother
        father = template.father
        suspected_father = template.suspected_father
        suspected_mother = template.suspected_mother
        sex = template.sex
        node = self.generate_node(father, mother, suspected_father,
                                  suspected_mother, sex = sex,
                                  twin = template)
        template.set_twin(node)
        return node
    
    @property
    def mapping(self):
        return self._mapping

class Node:
    """
    The node class represents a person in our genealogy.
    The mother and father properties are the true biological mother and father

    The suspected_mother and suspected_father are the individuals the
  
  'attacker' thinks are the mother and
    father. suspected_mother/father may be the true biological mother
    and father, or they may be other nodes in the geneology. They can
    be different due to errors in the records the attacker has, which
    can be caused by non-paternity events or adoption for example.
    """
    def __init__(self, node_generator, self_id, father_id = None,
                 mother_id = None, suspected_father_id = None,
                 suspected_mother_id = None, sex = None, twin_id = None):
        self.genome = None
        self._suspected_genome = None
        self._node_generator = node_generator
        self._twin_id = twin_id
        self._mother_id = mother_id
        if suspected_mother_id is None:
            self._suspected_mother_id = mother_id
        else:
            self._suspected_mother_id = suspected_mother_id
        if suspected_father_id is None:
            self._suspected_father_id = father_id
        else:
            self._suspected_father_id = suspected_father_id
        self._father_id = father_id
        self._id = self_id
        if isinstance(sex, Sex):
            self.sex = sex
        else:
            self.sex = choice(SEXES)
        self._children = []
        self._suspected_children = []
        self._resolve_parents()
        if self.mother is not None:
            assert self.mother.sex == Sex.Female
            self.mother._children.append(self._id)
        if self.father is not None:
            assert self.father.sex == Sex.Male
            self.father._children.append(self._id)
        self.set_suspected_mother(self.suspected_mother)
        self.set_suspected_father(self.suspected_father)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["mother"]
        del state["father"]
        del state["suspected_mother"]
        del state["suspected_father"]
        del state["twin"]
        return state
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        if "_suspected_genome" not in state:
            self._suspected_genome = None #We do this because there are old pickle files that don't have this property
        # Parents are resolved by PopulationUnpickler after all Node
        # objects are created

    def __str__(self):
        to_str = [STR_BASE.format(self._id, self._mother_id, self._father_id)]
        if self._twin_id is not None:
            to_str.append(", twin id = {}".format(self._twin_id))
        if self._suspected_mother_id != self._mother_id:
            to_str.append(", suspected mother id = {}".format(self._suspected_mother_id))
        if self._suspected_father_id != self._father_id:
            to_str.append(", suspected father id = {}".format(self._suspected_father_id))
        to_str.append(">")
        return "".join(to_str)


    def _resolve_parents(self):
        """
        This method resolves the respective parent ids in to objects.
        This is usually done in the constructor or after unpickling a node
        object.
        """
        mapping = self._node_generator._mapping

        if self._twin_id is not None:
            self.twin = mapping[self._twin_id]
        else:
            self.twin = None

        if self._mother_id is not None:
            self.mother = mapping[self._mother_id]
        else:
            self.mother = None
        if self._father_id is not None:
            self.father = mapping[self._father_id]
        else:
            self.father = None
            
        if self._suspected_father_id is not None:
            self.suspected_father = mapping[self._suspected_father_id]
        else:
            self.suspected_father = None
        if self._suspected_mother_id is not None:
            self.suspected_mother = mapping[self._suspected_mother_id]
        else:
            self.suspected_mother = None

    def set_twin(self, twin):
        self._twin_id = twin._id
        self.twin = twin

    # TODO: Turn suspected mother and father into @property methods
    def set_suspected_mother(self, suspected_mother):
        if (self.suspected_mother is not None
            and self._id in self.suspected_mother._suspected_children):
            self.suspected_mother._suspected_children.remove(self._id)
        self.suspected_mother = suspected_mother
        if suspected_mother is not None:
            assert suspected_mother.sex == Sex.Female
            self._suspected_mother_id = suspected_mother._id
            suspected_mother._suspected_children.append(self._id)
        else:
            self._suspected_mother_id = None

    def set_suspected_father(self, suspected_father):
        if (self.suspected_father is not None
            and self._id in self.suspected_father._suspected_children):
            self.suspected_father._suspected_children.remove(self._id)
        self.suspected_father = suspected_father
        if suspected_father is not None:
            assert suspected_father.sex == Sex.Male
            self._suspected_father_id = suspected_father._id
            suspected_father._suspected_children.append(self._id)
        else:
            self._suspected_father_id = None

    @property
    def suspected_genome(self):
        if self._suspected_genome is None:
            return self.genome
        else:
            return self._suspected_genome

    @suspected_genome.setter
    def suspected_genome(self, genome):
        self._suspected_genome = genome
            
    @property
    def mapping(self):
        """
        Returns a the dictionary mapping node id -> node object
        """
        return self._node_generator._mapping

    @property
    def children(self):
        """
        The true children of this node
        """
        return [self._node_generator._mapping[node_id]
                for node_id in self._children]

    @property
    def suspected_children(self):
        """
        The suspected children of this node
        """
        return [self._node_generator._mapping[node_id]
                for node_id in self._suspected_children]

    @property
    def node_generator(self):
        return self._node_generator
