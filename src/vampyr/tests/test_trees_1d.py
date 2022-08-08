from vampyr import BottomUp, Hilbert, Lebesgue, TopDown
from vampyr import vampyr1d as vp

r0 = [0.1]
r1 = [-0.1]

D = 1
k = 5
N = -1
n = 1
l = [2]
name = "func"
two_d = 2**D
kp1_d = (k + 1) ** D

world = vp.BoundingBox(scale=N)
root = vp.NodeIndex(scale=N)
idx = vp.NodeIndex(scale=n, translation=l)
mra = vp.MultiResolutionAnalysis(box=world, order=k)
tree = vp.FunctionTree(mra, name)


def test_FunctionTree():
    assert tree.norm() < 0.0
    assert tree.squaredNorm() < 0.0
    assert tree.nNodes() == 1
    assert tree.nEndNodes() == 1
    assert tree.nGenNodes() == 0
    assert tree.rootScale() == N
    assert tree.depth() == 1
    assert tree.name() == name
    assert tree.MRA() == mra


def test_FunctionTreeZero():
    tree.setZero()

    assert tree.norm() == 0.0
    assert tree.squaredNorm() == 0.0
    assert tree.integrate() == 0.0
    assert tree(r0) == 0.0
    assert tree(r1) == 0.0

    tree.clear()

    assert tree.norm() < 0.0
    assert tree.squaredNorm() < 0.0


def test_FunctionTreeSave():
    tree.setZero()
    path = tree.saveTree(filename=name)

    tree_2 = vp.FunctionTree(mra)
    tree_2.setName("func_2")
    tree_2.loadTree(filename=name)

    assert tree_2.name() == "func_2"
    assert tree_2.norm() == 0.0
    assert tree_2.squaredNorm() == 0.0
    assert tree_2.integrate() == 0.0
    assert tree_2.MRA() == tree.MRA()
    assert tree_2(r0) == 0.0
    assert tree_2(r1) == 0.0

    path.unlink()


def test_NodeIndex():
    parent = root.parent()
    child_0 = root.child(0)
    child_1 = root.child(1)

    assert root.scale() == N
    assert root.translation() == [0]
    assert parent.scale() == N - 1
    assert child_0.scale() == N + 1
    assert child_1.scale() == N + 1
    assert parent.translation() == [0]
    assert child_0.translation() == [0]
    assert child_1.translation() == [1]
    assert child_0 != root
    assert child_0 != parent
    assert child_0 != child_1
    assert child_0.parent() == root
    assert child_1.parent() == root
    assert parent.child(0) == root
    assert parent.child(1) != root

    parent.setScale(N)
    child_1.setTranslation([0])

    assert parent == root
    assert child_1 == child_0


def test_RootNode():
    tree.setZero()

    for i in range(tree.nRootNodes()):
        root_i = tree.fetchRootNode(i)
        assert root_i.hasCoefs()
        assert root_i.isAllocated()
        assert root_i.isRootNode()
        assert root_i.isEndNode()
        assert root_i.isLeafNode()
        assert not root_i.isBranchNode()
        assert not root_i.isGenNode()
        assert not root_i.hasParent()
        assert root_i.nChildren() == 0
        assert root_i.depth() == 0
        assert root_i.scale() == N
        assert root_i.nCoefs() == two_d * kp1_d
        assert root_i.integrate() == 0.0
        assert root_i.squaredNorm() == 0.0
        assert root_i.scalingNorm() == 0.0
        assert root_i.waveletNorm() == 0.0
        for t in range(two_d):
            assert root_i.componentNorm(t) == 0.0


def test_GenNode():
    tree.setZero()

    assert tree.nNodes() == 1
    assert tree.nEndNodes() == 1
    assert tree.nGenNodes() == 0

    node = tree.fetchNode(idx)

    assert tree.nNodes() == 1
    assert tree.nEndNodes() == 1
    assert tree.nGenNodes() == 4

    two_d = 2**D
    kp1_d = (k + 1) ** D
    assert node.hasCoefs()
    assert node.isAllocated()
    assert node.isLeafNode()
    assert node.isGenNode()
    assert node.hasParent()
    assert not node.isEndNode()
    assert not node.isRootNode()
    assert not node.isBranchNode()
    assert node.nChildren() == 0
    assert node.depth() == n - N
    assert node.scale() == n
    assert node.nCoefs() == kp1_d
    assert node.integrate() == 0.0
    assert node.squaredNorm() == 0.0
    assert node.scalingNorm() == 0.0
    assert node.waveletNorm() == 0.0
    for t in range(two_d):
        assert node.componentNorm(t) == 0.0

    tree.deleteGenerated()

    assert tree.nNodes() == 1
    assert tree.nEndNodes() == 1
    assert tree.nGenNodes() == 0


def test_HilbertIterator():
    tree.setZero()
    tree.fetchNode(idx)  # generate extra nodes

    it = vp.TreeIterator(tree, iterator=Hilbert)
    it.setMaxDepth(-1)
    it.setTraverse(TopDown)
    it.setReturnGenNodes(False)

    gen_count = 0
    node_count = 0
    while it.next():
        node = it.get()
        if node.isGenNode():
            gen_count += 1
        else:
            node_count += 1

    assert gen_count == 0
    assert node_count == tree.nNodes()


def test_LebesgueIterator():
    tree.setZero()
    tree.fetchNode(idx)  # generate extra nodes

    it = vp.TreeIterator(traverse=BottomUp)
    it.setMaxDepth(-1)
    it.setIterator(Lebesgue)
    it.setReturnGenNodes(True)

    gen_count = 0
    node_count = 0

    it.init(tree)
    while it.next():
        node = it.get()
        if node.isGenNode():
            gen_count += 1
        else:
            node_count += 1

    assert gen_count == tree.nGenNodes()
    assert node_count == tree.nNodes()
