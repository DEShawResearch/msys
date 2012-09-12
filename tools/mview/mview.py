#!/usr/bin/env desres-exec
#{
# desres-cleanenv --user-env \
# -m Python/2.7.1-06A/bin \
# -m numpy/1.5.1-29A/lib-python \
# -m msys/1.6.2/lib-python \
# -m sip/4.9-03A/lib-python \
# -m PyQt4/4.6-03A/lib-python \
# -- python $0 "$@"
#}

from PyQt4 import QtCore, QtGui
Qt = QtCore.Qt

import msys

class MsysTreeModel(QtCore.QAbstractItemModel):
    def __init__(self, mol, parent=None):
        super(MsysTreeModel, self).__init__(parent)
        self.rootItem = SystemTreeItem(None, mol)

    def columnCount(self, parent=QtCore.QModelIndex()):
        #print("columnCount")
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        #print("data")
        if not index.isValid():
            return QtCore.QVariant()

        if role != Qt.DisplayRole:
            return QtCore.QVariant()

        item = index.internalPointer()
        return item.data(index.column())

    def flags(self, index):
        #print("flags")
        if not index.isValid():
            return 0
        return QtCore.Qt.ItemFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)

    def headerData(self, section, orientation, role):
        #print("headerData")
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.rootItem.data(section)
        return QtCore.QVariant()

    def index(self, row, column, parent):
        #print("index")
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem is not None:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def parent(self, index):
        #print("parent")
        if not index.isValid():
            return QtCore.QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parent()

        if parentItem is self.rootItem:
            return QtCore.QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent):
        #print("rowCount")
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount()

class TreeItem(object):
    def __init__(self,parent,data):
        self.parentItem = parent
        self.itemData = data
        self.childItems = []

    def child(self,row):
        return self.childItems[row]

    def childCount(self):
        return len(self.childItems)

    def columnCount(self):
        return len(self.itemData)
 
    def data(self,column):
        return self.itemData[column]

    def parent(self):
        return self.parentItem

    def row(self):
        if self.parentItem is not None:
            return self.parentItem.childItems.index(self)
        return 0

def SystemTreeItem(parent, mol):
    tree = TreeItem(parent, [mol.name])
    tree.childItems = [ParticlesTreeItem(tree, mol), TablesTreeItem(tree, mol)]
    return tree

def ParticlesTreeItem(parent, mol):
    tree = TreeItem(parent, ["Particles"])
    for c in mol.chains:
        tree.childItems.append(ChainTreeItem(tree, c))
    return tree

def TablesTreeItem(parent, mol):
    tree = TreeItem(parent, ["Tables"])
    for table in mol.tables:
        tree.childItems.append(TableTreeItem(tree, table))
    return tree

def ChainTreeItem(parent, chain):
    tree = TreeItem(parent, ["Chain '%s'" % chain.name])
    for r in chain.residues:
        tree.childItems.append(ResidueTreeItem(tree, r))
    return tree

def ResidueTreeItem(parent, residue):
    tree = TreeItem(parent, ["%4s %d" % (residue.name, residue.resid)])
    for a in residue.atoms:
        tree.childItems.append(AtomTreeItem(tree, a))
    return tree

def AtomTreeItem(parent, atom):
    tree = TreeItem(parent, ["%d %s" % (atom.id, atom.name)])
    return tree

def TableTreeItem(parent, table):
    tree = TreeItem(parent, [table.name])
    return tree

if __name__=="__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    for path in sys.argv[1:]:
        print "Loading", path
        mol=msys.Load(path)
        view = QtGui.QTreeView()
        print "Processing..."
        model = MsysTreeModel(mol)
        print "done"
        view.setModel(model)
        view.setWindowTitle("MView")
        view.show()

    sys.exit(app.exec_())




