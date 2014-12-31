import re

class inode(object):
	def __init__(self, name, parent, data):
		self.__name = str(name)
		self.__children = []
		self.__parent = None
		self.data = data
		if isinstance(parent, inode) :
			parent.add(self)
	def __getitem__(self, index) :
		return self.__children[index]
	def __contains__(self, node) :
		return node in self.__children
	def __len__(self) :
		return len( self.__children )
	def __iter__(self) :
		return self.dfs()
	def is_root(self) :
		return self.__parent is None
	def parent(self) :
		return self.__parent;
	def name(self) :
		return self.__name;
	def path(self, node = None) :
		path = [ self ]
		while not path[0].is_root() and path[0] != node :
			path = [ path[0].__parent ] + path
		return path
	def unlink(self) :
		if not self.is_root() :
			self.__parent.__children.remove(self)
			self.__parent = None
		return self
	def __link(self, node) :
		self.__children.append(node)
		node.__parent = self
		return node
	def add(self, node) :
		if node not in self.path() :
			self.__link(node.unlink())
		return self
	def rem(self, node) :
		if node in self :
			node.unlink()
		return node
	def clear(self) :
		stack, queue = [self], []
		while stack :
			node = stack.pop(-1)
			if not node.__children:
				node.__parent = None
				queue.append(node)
				continue
			stack.append(node)
			stack.extend([child for child in node.__children])
			node.__children = []
		return queue
	def fetch(self, name) :
		for node in self.__children :
			if node.__name == name :
				return node
		return None 
	def children(self, names = None) :
		if names is None :
			return self.__children
		return [node for node in self.__children if node.__name in names]
	def find(self, pattern) :
		if pattern is None :
			return self.__children
		try :
			regex = re.compile(pattern, re.UNICODE)
			return [node for node in self.__children
				if regex.search(node.__name)]
		except :
			return []
	def descends_from(self, node) :
		return node is not self and node in self.path()
	## Locate the last common ancestor
	def LCA_with(self, node) :
		## Get the lineages
		lin0, lin1 = self.path(), node.path()
		## Find the longest common suffix
		inx = [i for i in range(min(len(lin0), len(lin1))) if lin0[i] is lin1[i]]
		return lin0[max(inx)] if inx else None
	def dfs(self) :
		stack = [(self, 0, False)]
		while stack :
			node, level, visited = stack.pop(-1)
			if visited or not node.__children:
				yield node, level
				continue
			stack.append((node, level, True))
			stack.extend([(child, level+1, False) for child in node.__children])
	def dfs_tree(self) :
		queue = [ (self, 0) ]
		while queue :
			node, level = queue.pop(0)
			yield node, level
			queue = [(child,level+1) for child in node.__children] + queue
	def bfs(self) :
		queue = [ (self, 0) ]
		while queue :
			node, level = queue.pop(0)
			yield node, level
			queue.extend([(child,level+1) for child in node.__children])
	def dfs_recur(self, level=0) :
		for node in self.__children :
			for x in node.dfs(level+1) :
				yield x
		yield (self, level)
	def uri(self) :
		return "/".join([node.__name for node in self.path()])
	def __repr__(self) :
		return self.uri()
	def rename(self, new_name) :
		old_name = self.__name
		self.__name = new_name
		return old_name

"""Usage example"""
"""
tree = inode("*", None, {}).add(
	inode("a", None, {}).add(
		inode("a", None, {})).add(
		inode("b", None, {}).add(
			inode("a", None, {})))).add(
	inode("b", None, {}))
[x for x in tree]
[x for x in tree.dfs_tree()]

tree[0][0].path()
tree[0][0].LCA_with(tree[0][1][0])
tree[0][0].LCA_with(tree[1])

[x for x in tree]
ch = tree.clear()

[x for x in tree]
ch
ch[2].path()
"""
