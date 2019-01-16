from sex import Sex

class Generation:
    def __init__(self, members = None):
        self.members = list(members)

    @property
    def men(self):
        return (person for person in self.members if person.sex == Sex.Male)

    @property
    def women(self):
        return (person for person in self.members if person.sex == Sex.Female)
    
    @property
    def size(self):
        return len(self.members)
