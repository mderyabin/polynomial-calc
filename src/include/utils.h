#ifndef __UTILS_H__
#define __UTILS_H__

#include <fstream>
#include <string>

#include <cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

namespace polycalc {

enum SER_Archive_Type {BIN, JSON};

template<typename T>
bool SerializeObject(T &obj, std::string filename, const SER_Archive_Type TYPE = BIN) {
    std::fstream fs;

    switch (TYPE) {
    case JSON: {
        fs.open(filename, std::ios::out);
        if (fs.is_open()) { 
            { // to finish json inside brakets... 
                cereal::JSONOutputArchive jarchive(fs);
                jarchive(obj);
            }
            fs.close();
            return true;
        } else {
            return false;
        }
        break;
    }
    case BIN: {
        fs.open(filename, std::ios::out | std::ios::binary);
        if (fs.is_open()) {
            {
                cereal::BinaryOutputArchive barchive(fs);
                barchive(obj);
            }
            fs.close();
            return true;
        } else {
            return false;
        }
        break;
    }
    default: 
        break;
    }

    return false;
}

template<typename T>
bool DeserializeObject(T &obj, std::string filename, const SER_Archive_Type TYPE = BIN) {
    std::fstream fs;

    switch (TYPE) {
    case JSON: {
        fs.open(filename, std::ios::in);
        if (fs.is_open()) {
            cereal::JSONInputArchive jarchive(fs);
            jarchive(obj);
            fs.close();
            return true;
        } else {
            return false;
        }
        break;
    }
    case BIN: {
        fs.open(filename, std::ios::in | std::ios::binary);
        if (fs.is_open()) {
            cereal::BinaryInputArchive barchive(fs);
            barchive(obj);
            fs.close();
            return true;
        } else {
            return false;
        }
        break;
    }
    default: 
        break;
    }

    return false;
}

}


#endif /* __UTILS_H__ */