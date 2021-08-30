import numpy as np
import vampyr as vp
import pytest

r0 = [0.1, 0.1, 0.1]
r1 = [-0.1, -0.1, -0.1]

D = 3
k = 5
N = -1
n = 1
l = [2, 2, 2]
name = "func"
two_d = 2**D
kp1_d = (k + 1)**D

legendre = vp.LegendreBasis(order=k)
world = vp.D3.BoundingBox(scale=N)
root = vp.D3.NodeIndex(scale=N)
idx = vp.D3.NodeIndex(scale=n, translation=l)
mra = vp.D3.MultiResolutionAnalysis(box=world, basis=legendre)
tree = vp.D3.FunctionTree(mra)
tree.setName(name)

def test_FunctionTree():
    assert tree.getSquareNorm() < 0.0
    assert tree.getNNodes() == 1
    assert tree.getNEndNodes() == 1
    assert tree.getNGenNodes() == 0
    assert tree.getRootScale() == N
    assert tree.getDepth() == 1
    assert tree.getName() == name
    assert tree.getMRA() == mra

def test_FunctionTreeZero():
    tree.setZero()

    assert tree.getSquareNorm() == 0.0
    assert tree.integrate() == 0.0
    assert tree(r0) == 0.0
    assert tree(r1) == 0.0

    tree.clear()

    assert tree.getSquareNorm() < 0.0

def test_FunctionTreeSave():
    tree.setZero()
    tree.saveTree(filename=name)

    tree_2 = vp.D3.FunctionTree(mra)
    tree_2.setName("func_2")
    tree_2.loadTree(filename=name)

    assert tree_2.getName() == "func_2"
    assert tree_2.getSquareNorm() == 0.0
    assert tree_2.integrate() == 0.0
    assert tree_2.getMRA() == tree.getMRA()
    assert tree_2(r0) == 0.0
    assert tree_2(r1) == 0.0

def test_NodeIndex():
    parent = root.parent()
    child_0 = root.child(0)
    child_1 = root.child(1)

    assert root.getScale() == N
    assert root.getTranslation() == [0, 0, 0]
    assert parent.getScale() == N - 1
    assert child_0.getScale() == N + 1
    assert child_1.getScale() == N + 1
    assert parent.getTranslation() == [0, 0, 0]
    assert child_0.getTranslation() == [0, 0, 0]
    assert child_1.getTranslation() == [1, 0, 0]
    assert child_0 != root
    assert child_0 != parent
    assert child_0 != child_1
    assert child_0.parent() == root
    assert child_1.parent() == root
    assert parent.child(0) == root
    assert parent.child(1) != root

    parent.setScale(N)
    child_1.setTranslation([0, 0, 0])

    assert parent == root
    assert child_1 == child_0

def test_RootNode():
    tree.setZero()

    for i in range(tree.getNRootNodes()):
        root_i = tree.getRootNode(i)
        assert root_i.hasCoefs()
        assert root_i.isAllocated()
        assert root_i.isRootNode()
        assert root_i.isEndNode()
        assert root_i.isLeafNode()
        assert not root_i.isBranchNode()
        assert not root_i.isGenNode()
        assert not root_i.hasParent()
        assert root_i.getNChildren() == 0
        assert root_i.getDepth() == 0
        assert root_i.getScale() == N
        assert root_i.getNCoefs() == two_d * kp1_d
        assert root_i.integrate() == 0.0
        assert root_i.getSquareNorm() == 0.0
        assert root_i.getScalingNorm() == 0.0
        assert root_i.getWaveletNorm() == 0.0
        for t in range(two_d):
            assert root_i.getComponentNorm(t) == 0.0

def test_GenNode():
    tree.setZero()

    assert tree.getNNodes() == 1
    assert tree.getNEndNodes() == 1
    assert tree.getNGenNodes() == 0

    node = tree.getNode(idx)

    assert tree.getNNodes() == 1
    assert tree.getNEndNodes() == 1
    assert tree.getNGenNodes() == 16

    two_d = 2**D
    kp1_d = (k + 1)**D
    assert node.hasCoefs()
    assert node.isAllocated()
    assert node.isLeafNode()
    assert node.isGenNode()
    assert node.hasParent()
    assert not node.isEndNode()
    assert not node.isRootNode()
    assert not node.isBranchNode()
    assert node.getNChildren() == 0
    assert node.getDepth() == n - N
    assert node.getScale() == n
    assert node.getNCoefs() == kp1_d
    assert node.integrate() == 0.0
    assert node.getSquareNorm() == 0.0
    assert node.getScalingNorm() == 0.0
    assert node.getWaveletNorm() == 0.0
    for t in range(two_d):
        assert node.getComponentNorm(t) == 0.0

    tree.deleteGenerated()

    assert tree.getNNodes() == 1
    assert tree.getNEndNodes() == 1
    assert tree.getNGenNodes() == 0

def test_HilbertIterator():
    tree.setZero()
    tree.getNode(idx) # generate extra nodes

    it = vp.D3.HilbertIterator(tree)
    it.setMaxDepth(-1)
    it.setReturnGenNodes(False)

    gen_count = 0
    node_count = 0
    while (it.next()):
        node = it.getNode()
        if node.isGenNode():
            gen_count += 1
        else:
            node_count += 1

    assert gen_count == 0
    assert node_count == tree.getNNodes()

def test_LebesgueIterator():
    tree.setZero()
    tree.getNode(idx) # generate extra nodes

    it = vp.D3.LebesgueIterator(tree)
    it.setMaxDepth(-1)
    it.setReturnGenNodes(True)

    gen_count = 0
    node_count = 0
    while (it.next()):
        node = it.getNode()
        if node.isGenNode():
            gen_count += 1
        else:
            node_count += 1

    assert gen_count == tree.getNGenNodes()
    assert node_count == tree.getNNodes()
