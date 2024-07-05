#include <iostream>
#include <string>
using namespace std;

void Swap(string &string1, string &string2);

int main()
{   
    string string1 = "Kool - Aid";
    string string2 = "Water";
    string temp;

    Swap(string1, string2);
    

    cout << "X: " << string1 << '\n';
    cout << "Y: " << string2 << '\n';
    
    return 0;


}

void Swap(string &string1, string &string2){
    string temp;
    
    temp = string1;
    string1 = string2;
    string2 = temp;
    
}