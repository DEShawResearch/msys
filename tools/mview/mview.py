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

class SystemTreeItem(TreeItem):
    def __init__(self, parent, mol):
        TreeItem.__init__(self, parent, [mol.name, ''])
        self.childItems = [
                ParticlesTreeItem(self, mol), 
                TablesTreeItem(self, mol),
                ]

class ParticlesTreeItem(TreeItem):
    def __init__(self, parent, mol):
        TreeItem.__init__(self, parent, ["Particles", str(mol.natoms)])
        for c in mol.chains:
            self.childItems.append(ChainTreeItem(self, c))

class TablesTreeItem(TreeItem):
    def __init__(self, parent, mol):
        name="Tables"
        size=str(len(mol.tables))
        TreeItem.__init__(self, parent, [name,size])
        for table in mol.tables:
            self.childItems.append(TableTreeItem(self, table))

class ChainTreeItem(TreeItem):
    def __init__(self, parent, chain):
        name="Chain '%s'" % chain.name
        size=str(chain.nresidues)
        TreeItem.__init__(self, parent, [name, size])
        for r in chain.residues:
            self.childItems.append(ResidueTreeItem(self, r))

class ResidueTreeItem(TreeItem):
    def __init__(self, parent, residue):
        name="%4s %d" % (residue.name, residue.resid)
        size=str(residue.natoms)
        TreeItem.__init__(self, parent, [name, size])
        for a in residue.atoms:
            self.childItems.append(AtomTreeItem(self, a))

class AtomTreeItem(TreeItem):
    def __init__(self, parent, atom):
        name="%d %s" % (atom.id, atom.name)
        size=''
        TreeItem.__init__(self, parent, [name, size])

class TableTreeItem(TreeItem):
    def __init__(self, parent, table):
        name=table.name
        size="%d params, %d terms" % (table.params.nparams, table.nterms)
        TreeItem.__init__(self, parent, [name,size])

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




