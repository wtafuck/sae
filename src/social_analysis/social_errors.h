#pragma once
#include <iostream>
#include <string>
class Exception
{
    std::string info;
public:
    Exception(std::string _info):info(_info){}
    void print(){
        std::cerr << "Error: "<< info <<std::endl;
    }
};