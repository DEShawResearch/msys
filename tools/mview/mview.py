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
    def __init__(self, parent=None):
        super(MsysTreeModel, self).__init__(parent)
        self.rootItem = SystemListItem(None)

    def addSystem(self, mol):
        self.rootItem.addChild(mol)

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

    def row(self):
        if self.parentItem is not None:
            return self.parentItem.childItems.index(self)
        return 0

    def columnCount(self):
        return len(self.itemData)
 
    def data(self,column):
        return self.itemData[column]

    def parent(self):
        return self.parentItem

class SystemListItem(TreeItem):
    def __init__(self, parent):
        TreeItem.__init__(self, parent, ["Systems", ""])

    def addChild(self, mol):
        self.childItems.append(SystemTreeItem(self, mol))

class SystemTreeItem(TreeItem):
    def __init__(self, parent, mol):
        TreeItem.__init__(self, parent, [mol.name, "%d atoms" % mol.natoms])
        self.childItems = [
                ParticlesTreeItem(self, mol), 
                TablesTreeItem(self, mol),
                ]

class ParticlesTreeItem(TreeItem):
    def __init__(self, parent, mol):
        TreeItem.__init__(self, parent, ["Structure", "%d chains" % mol.nchains])
        for c in mol.chains:
            self.childItems.append(ChainTreeItem(self, c))

class TablesTreeItem(TreeItem):
    def __init__(self, parent, mol):
        name="Forcefield"
        size="%d tables" % len(mol.tables)
        TreeItem.__init__(self, parent, [name,size])
        for table in mol.tables:
            self.childItems.append(TableTreeItem(self, table))

class ChainTreeItem(TreeItem):
    def __init__(self, parent, chain):
        name="Chain '%s'" % chain.name
        size="%d residues" % chain.nresidues
        TreeItem.__init__(self, parent, [name, size])
        for r in chain.residues:
            self.childItems.append(ResidueTreeItem(self, r))

class ResidueTreeItem(TreeItem):
    def __init__(self, parent, residue):
        name="%4s %d" % (residue.name, residue.resid)
        size="%d atoms" % residue.natoms
        TreeItem.__init__(self, parent, [name, size])
        self.residue = residue
        self.childItems = None

    def build(self):
        ''' lazy building of atoms '''
        if self.childItems is None:
            self.childItems = [AtomTreeItem(self,a) for a in self.residue.atoms]

    def child(self,row):
        self.build()
        return self.childItems[row]

    def childCount(self):
        return self.residue.natoms


class AtomTreeItem(TreeItem):
    def __init__(self, parent, atom):
        name="%d %s" % (atom.id, atom.name)
        size=''
        TreeItem.__init__(self, parent, [name, size])

class TableTreeItem(TreeItem):
    def __init__(self, parent, table):
        name=table.name
        size=table.category
        #size="%d params, %d terms" % (table.params.nparams, table.nterms)
        TreeItem.__init__(self, parent, [name,size])
        self.childItems = [
                ParamsItem(self, table.params),
                TermsItem(self, table),
                ]

class TermsItem(TreeItem):
    def __init__(self, parent, table):
        name="Terms"
        size=table.nterms
        TreeItem.__init__(self, parent, [name,size])
        self.table = table
        self.childItems = None

    def build(self):
        ''' lazy building of terms '''
        if self.childItems is None:
            self.childItems = [TermTreeItem(self,t) for t in self.table.terms]

    def child(self,row):
        self.build()
        return self.childItems[row]

    def childCount(self):
        return self.table.nterms

class TermTreeItem(TreeItem):
    def __init__(self, parent, term):
        name="%d %s" % (term.id, [a.id for a in term.atoms])
        size=[]
        for p in term.table.params.props:
            size.append('%s=%s' % (p, term[p]))
        size=' '.join(size)
        TreeItem.__init__(self, parent, [name,size])

class ParamsItem(TreeItem):
    def __init__(self, parent, params):
        name="Params"
        size=params.nparams
        TreeItem.__init__(self, parent, [name,size])
        self.params = params
        self.childItems = None

    def build(self):
        ''' lazy building of params '''
        if self.childItems is None:
            self.childItems = [ParamItem(self,p) for p in self.params.params]

    def child(self,row):
        self.build()
        return self.childItems[row]

    def childCount(self):
        return self.params.nparams

class ParamItem(TreeItem):
    def __init__(self, parent, param): 
       name=param.id
       size=[]
       for p in param.table.props:
           size.append('%s=%s' % (p,param[p]))
       size=' '.join(size)
       TreeItem.__init__(self, parent, [name,size])

if __name__=="__main__":
    import sys
    model = MsysTreeModel()
    for path in sys.argv[1:]:
        print "Loading", path
        mol=msys.Load(path)
        print "Processing..."
        model.addSystem(mol)
        print "done"

    app = QtGui.QApplication(sys.argv)
    view = QtGui.QTreeView()
    view.setMinimumSize(500,500)
    view.setModel(model)
    view.resizeColumnToContents(0)
    win=QtGui.QMainWindow()
    win.setWindowTitle("MView")
    win.setCentralWidget(view)
    win.show()
    sys.exit(app.exec_())




