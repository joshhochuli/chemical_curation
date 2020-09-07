import unittest

import sys
sys.path.append("../src")
from chemical_curation import modification_graph

class TestModificationGraph(unittest.TestCase):

    def test_modification_base_init(self):
        mod = modification_graph.Modification(text = "modification made")
        print('\n', mod)

    def test_modification_full_init(self):
        parent_mod = modification_graph.Modification(text = "parent")
        child_mod = modification_graph.Modification(text = "child")
        mod = modification_graph.Modification(text = "main", parent = parent_mod, child = child_mod)
        self.assertEqual(mod.parents[0], parent_mod)
        self.assertEqual(mod.child,child_mod)

    def test_replace_parent(self):
        parent_mod_1 = modification_graph.Modification(text = "parent1")
        parent_mod_2 = modification_graph.Modification(text = "parent2")
        mod = modification_graph.Modification(text = "main", parent = parent_mod_2)
        mod.set_parent(parent_mod_2, append = False)
        self.assertEqual(mod.parents[0],parent_mod_2)

    def test_append_parent(self):
        parent_mod_1 = modification_graph.Modification(text = "parent1")
        parent_mod_2 = modification_graph.Modification(text = "parent2")
        mod = modification_graph.Modification(text = "main", parent = parent_mod_1)
        mod.set_parent(parent_mod_2, append = True)
        self.assertEqual(len(mod.parents), 2)
        self.assertEqual(mod.parents[0],parent_mod_1)
        self.assertEqual(mod.parents[1],parent_mod_2)

    def test_init_graph(self):
        mod = modification_graph.Modification(text = "first")
        graph = modification_graph.Modification_Graph(head = mod)

    def test_add_graph(self):
        mod1 = modification_graph.Modification(text = "first")
        mod2 = modification_graph.Modification(text = "second")
        graph = modification_graph.Modification_Graph(head = mod1)
        graph.add_modification(mod2)
        self.assertEqual(graph.head.parents[0], mod1)
        self.assertEqual(len(graph.head.parents), 1)
        self.assertEqual(graph.head.parents[0].parents, None)
        self.assertEqual(graph.head.child, None)
        self.assertEqual(graph.head.parents[0].child, mod2)

    def test_print_graph(self):
        mod1 = modification_graph.Modification(text = "first")
        mod2 = modification_graph.Modification(text = "second")
        mod3 = modification_graph.Modification(text = "third")
        graph = modification_graph.Modification_Graph(head = mod1)
        graph.add_modification(mod2)
        graph.add_modification(mod3)
        graph.print()

    def test_graph_len(self):
        mod1 = modification_graph.Modification(text = "first")
        mod2 = modification_graph.Modification(text = "second")
        mod3 = modification_graph.Modification(text = "third")
        graph = modification_graph.Modification_Graph(head = mod1)
        self.assertEqual(len(graph), 1)
        graph.add_modification(mod2)
        self.assertEqual(len(graph), 2)
        graph.add_modification(mod3)
        self.assertEqual(len(graph), 3)

    def test_single_merge(self):
        mod1a = modification_graph.Modification(text = "1a")
        mod1b = modification_graph.Modification(text = "1b")
        graph_1 = modification_graph.Modification_Graph()
        graph_1.add_modification(mod1a)
        graph_1.add_modification(mod1b)
        self.assertEqual(len(graph_1), 2)
        mod2a = modification_graph.Modification(text = "2a")
        mod2b = modification_graph.Modification(text = "2b")
        graph_2 = modification_graph.Modification_Graph()
        graph_2.add_modification(mod2a)
        graph_2.add_modification(mod2b)
        self.assertEqual(len(graph_2), 2)

        graph_1_len = len(graph_1)
        graph_2_len = len(graph_2)
        graph_1.merge(graph_2)

        nodelist, count = graph_1._depth_first_traversal(graph_1.get_head(), nodelist = [], count = 0)
        print(nodelist)
        #self.assertEqual(len(graph_1), graph_1_len + graph_2_len + 1)

    def test_dual_merge(self):
        mod1a = modification_graph.Modification(text = "1a")
        mod1b = modification_graph.Modification(text = "1b")
        graph_1 = modification_graph.Modification_Graph()
        graph_1.add_modification(mod1a)
        graph_1.add_modification(mod1b)
        self.assertEqual(len(graph_1), 2)

        mod2a = modification_graph.Modification(text = "2a")
        mod2b = modification_graph.Modification(text = "2b")
        graph_2 = modification_graph.Modification_Graph()
        graph_2.add_modification(mod2a)
        graph_2.add_modification(mod2b)
        self.assertEqual(len(graph_2), 2)

        mod3a = modification_graph.Modification(text = "3a")
        mod3b = modification_graph.Modification(text = "3b")
        mod3c = modification_graph.Modification(text = "3c")
        mod3d = modification_graph.Modification(text = "3d")
        graph_3 = modification_graph.Modification_Graph()
        graph_3.add_modification(mod3a)
        graph_3.add_modification(mod3b)
        graph_3.add_modification(mod3c)
        graph_3.add_modification(mod3d)


        graph_1_len = len(graph_1)
        graph_2_len = len(graph_2)
        graph_3_len = len(graph_3)
        graph_1.merge(graph_2)
        graph_1.merge(graph_3)

        nodelist, count = graph_1._depth_first_traversal(graph_1.get_head(), nodelist = [], count = 0)
        print(nodelist)
        self.assertEqual(len(graph_1), graph_1_len + graph_2_len + graph_3_len + 2)



if __name__ == '__main__':
    unittest.main()
