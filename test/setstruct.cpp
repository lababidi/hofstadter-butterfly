#include <iostream>
#include <set>
using namespace std;
struct foo 
{
    int key;
};

bool operator<(const foo& lhs, const foo& rhs)
{
    return lhs.key < rhs.key;
}

set<foo> bar;

int main()
{
    foo test, second;
    test.key = 0;
    second.key = 0;
    bar.insert(test);
    bar.insert(second);
    for(set<foo>::iterator iter=bar.begin(); iter!=bar.end();++iter) {
        cout<<(*iter).key<<" ";
    }
}
