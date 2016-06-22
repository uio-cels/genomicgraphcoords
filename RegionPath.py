class RegionPath(object):
    def __init__(self, id, metadata):
        self.id = id
        self.metadata = metadata
        self.linear_references = metadata

    def get_length(self):
        return max([lr.end-lr.start for lr in self.linear_references.values()])

    def contains(self, lin_ref):
        return any([my_ref.contains(lin_ref) for my_ref
                    in self.linear_references.values()])

    def intersects(self, lin_ref):
        return any([my_ref.intersects(lin_ref) for my_ref
                    in self.linear_references.values()])

    def length(self):
        return max([lr.length() for lr in self.linear_references.values()])

    def is_empty(self):
        if not self.linear_references:
            return True
        if (all([lr.start == lr.end for lr
                 in self.linear_references.values()])):
            return True
        return False

    def __str__(self):
        return " RegionPaths: " + ";".join([str(lin_ref) for lin_ref
                                            in self.linear_references.items()])

    def __repr__(self):
        return self.__str__()
