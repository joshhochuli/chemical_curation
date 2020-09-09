from datetime import datetime

class Modification_Graph:

    def __init__(self, head = None):

        self.head = head

    def set_head(self, head):

        self.head = head

    def get_head(self):

        return self.head

    def add_modification(self, mod = None, text = None):
        if mod != None and text != None:
            raise Exception("add_modification takes one of 'mod' and 'text', not both")

        if text:
            mod = Modification(text = text)

        if self.head:
            mod.set_parent(self.head)
            self.head.set_child(mod)

        self.head = mod

    def merge(self, other):

        merge_mod = Modification(text = "merge")
        merge_mod.set_parent(self.head)
        merge_mod.set_parent(other.get_head())
        other.get_head().set_child(merge_mod)
        self.head.set_child(merge_mod)
        self.head = merge_mod

    def print(self):

        curr_mod = self.head
        print('\n', curr_mod)
        while curr_mod.parents != None:
            curr_mod = curr_mod.parents[0]
            print(curr_mod)

    #returns True if any modification has `text` as a substring
    def has_modification(self, text):

        nodelist, count = self._depth_first_traversal(node = self.head, nodelist = [], count = 0)
        print(nodelist)
        found = False
        for node in nodelist:
            if type(node) == Modification:
                if text in node.text:
                    found = True
                    break

        return found

    def _depth_first_traversal(self, node, nodelist, count):

        nodelist.append(node)
        count += 1
        print("VISITING: ", node, count)

        if node.parents == None:
            return nodelist, count
        else:
            for parent in node.parents:
                if len(node.parents) > 1:
                    nodelist.append("BRANCH")
                ret_nodelist, ret_count = self._depth_first_traversal(parent, nodelist, count)
                #print("EXISTING: ",  nodelist)
                #print("FROM CALL: ", ret_val)
                nodelist = ret_nodelist
                count = ret_count

            return nodelist, count


    def __len__(self):

        if self.head == None:
            return 0
        else:
            '''
            count = 1
            curr_mod = self.head
            while curr_mod.parents != None:
                count += 1
                curr_mod = curr_mod.parents[0]

            return count
            '''
            nodelist, count = self._depth_first_traversal(self.head, [], 0)
            print(nodelist)
            return count

    def __repr__(self):
            nodelist, count = self._depth_first_traversal(self.head, [], 0)
            return "\n".join([str(node) for node in nodelist])




class Modification:
    ''' Stores a single modification to a chemical structure or associated data as a simple string.
    Acts as a node in a Modification_Graph
    Can have multiple parents (for merging structures) but only one child
    Examples: "Change =O to -OH", "'Na+' detected as salt and removed"

    '''

    def __init__(self, text, parent = None, child = None):
        self.text = text
        self.timestamp = datetime.now()
        if parent:
            self.parents = [parent]
        else:
            self.parents = None
        if child:
            self.child = child
        else:
            self.child = None

    def set_child(self, child):
        self.child = child

    def set_parent(self, parent, append = True):
        if not self.parents:
            self.parents = []
        if append:
            self.parents.append(parent)
        else:
            self.parents = [parent]

    def __repr__(self):

        return f"{self.timestamp}: {self.text}"

    def __str__(self):

        return self.__repr__()


