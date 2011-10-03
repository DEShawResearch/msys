#include "selection.hxx"
#include <string.h>

using namespace desres::msys;
using namespace desres::msys::atomsel;

Selection::Selection(Id size, const IdList& ids)
: flags(size) {
    for (unsigned i=0, n=ids.size(); i<n; i++) flags[i]=true;
}

void Selection::clear() {
    if (!flags.size()) return;
    memset(&flags[0], 0, flags.size()*sizeof(flags[0]));
}

Id Selection::count() const {
    Id sum=0;
    for (flag_list::const_iterator i=flags.begin(); i!=flags.end(); ++i) {
        sum += *i;
    }
    return sum;
}

Selection& Selection::intersect( const Selection& other ) {
    flag_list::iterator i=flags.begin();
    flag_list::iterator e=flags.end();
    flag_list::const_iterator j=other.flags.begin();
    for (; i!=e; ++i, ++j) {
        *i &= *j;
    }
    return *this;
}

Selection& Selection::add( const Selection& other ) {
    flag_list::iterator i=flags.begin();
    flag_list::iterator e=flags.end();
    flag_list::const_iterator j=other.flags.begin();
    for (; i!=e; ++i, ++j) {
        *i |= *j;
    }
    return *this;
}

Selection& Selection::subtract( const Selection& other ) {
    flag_list::iterator i=flags.begin();
    flag_list::iterator e=flags.end();
    flag_list::const_iterator j=other.flags.begin();
    for (; i!=e; ++i, ++j) {
        *i &= !(*j);
    }
    return *this;
}

IdList Selection::ids() const {
    IdList tmp;
    tmp.reserve(size());
    for (Id i=0, n=flags.size(); i!=n; i++) if (flags[i]) tmp.push_back(i);
    return tmp;
}
