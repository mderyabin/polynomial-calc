/**
 * @file utils.h
 * @author Msxim Deryabin (maxim.deryabin@gmail.com)
 * @brief Common utilities.
 * @version 0.1
 * @date 2023-05-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef __UTILS_H__
#define __UTILS_H__

#include <fstream>
#include <string>

#include <cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

namespace polycalc {

// Serialization type for cereal to dynamically choose the archive type
enum SER_Archive_Type {BIN, JSON};

/**
 * @brief Serialize general object of any type.
 * 
 * @tparam T serializable object (class or native type)
 * @param obj object to serialize
 * @param filename name of file with relative path
 * @param TYPE serialization type
 * @return true if serialized successfully
 * @return false file access error 
 */
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

/**
 * @brief Deserialize general object of any type.
 * 
 * @tparam T serializable object (class or native type)
 * @param obj object to deserialize
 * @param filename name of file with relative path
 * @param TYPE serialization type
 * @return true if deserialized successfully
 * @return false file access error 
 */
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