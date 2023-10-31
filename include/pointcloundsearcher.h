//
// Created by Zachary on 2023/8/28.
//

#ifndef PLYFASTLOADER_POINTCLOUNDSEARCHER_H
#define PLYFASTLOADER_POINTCLOUNDSEARCHER_H
#include "plyio.h"
#include <experimental/filesystem>
#include <memory>
#include <utility>
#include <json/json.h>
#include <unordered_map>
#include <vector>
#include <chrono>
#include "pct_entity.h"
#include <eigen3/Eigen/Dense>
#include <cmath>
namespace fs = std::experimental::filesystem::v1;
namespace PCT {
    class PointCloudSearcher {
    public:
        /**
         * use for create only
         * @param cachePath
         * @param depth
         * @param boundingBoxMin
         * @param boundingBoxMax
         */
        PointCloudSearcher(const fs::path& cachePath, int depth, const PCT::BoundingBox& boundingBox, bool override = false) {
            std::cout << "[" << currentTime() << "] create searcher for write start." << std::endl;

            if (cachePath.empty()) {
                throw std::runtime_error("[PCT]Cache path cannot be empty!");
            }
            if (fs::exists(cachePath)) {
                if (!fs::is_directory(cachePath)) {
                    throw std::runtime_error("[PCT]Invalid cache path:" + cachePath.string());
                }
            } else {
                fs::create_directories(cachePath);
            }

            this->cachePath = cachePath;
            this->pathIndexInfo = cachePath/"index.json";

            if (fs::exists(pathIndexInfo)) {
                if (override) {
                    fs::remove_all(cachePath);
                    fs::create_directories(cachePath);
                } else {
                    throw std::runtime_error("[PCT]Index info exists!");
                }
            }

            this->depth = depth;
            this->boundingBox = boundingBox;

            initIndexInfo();
            std::cout << "[" << currentTime() << "] create searcher for write end." << std::endl;
        }

        /**
         * use for read only
         * @param cachePath
         */
        explicit PointCloudSearcher(const fs::path& cachePath) {
            std::cout << "[" << currentTime() << "] create searcher for read start." << std::endl;
            if (cachePath.empty()) {
                throw std::runtime_error("[PCT]Cache path cannot be empty!");
            }
            if (fs::exists(cachePath)) {
                if (!fs::is_directory(cachePath)) {
                    throw std::runtime_error("[PCT]Invalid cache path:" + cachePath.string());
                }
            } else {
                throw std::runtime_error("[PCT]Cache path not exists!");
            }

            this->cachePath = cachePath;
            this->pathIndexInfo = cachePath/"index.json";

            readIndexInfo();
            std::cout << "[" << currentTime() << "] create searcher for read end." << std::endl;
        }

        ~PointCloudSearcher() = default;

        void addVertexData(PCT::PointCloudData& data) {
            if (data.getSpec().empty()) {
                throw std::runtime_error("[PCT]Spec in vertex data must be not empty!");
            }

            std::cout << "[" << currentTime() << "] add vertex data start." << std::endl;

            std::cout << "[" << currentTime() << "] start classify data with octree." << std::endl;

            std::unordered_map<std::string,std::unique_ptr<std::vector<size_t>>> octreeDataIndexMap;

            fastClassify(data, indexInfo["octree"], octreeDataIndexMap);

            std::cout << "[" << currentTime() << "] classify data done, begin write data and count points." << std::endl;

            size_t count = 0;
            for (auto& kv : octreeDataIndexMap) {
                count += kv.second->size();
            }
            std::cout << "[" << currentTime() << "] indexes count in map is:" << count << std::endl;

            for (auto& kv : octreeDataIndexMap) {
                writeSubPly(data, kv.first, kv.second);
                //统计点数
                Json::Value* octree = &indexInfo["octree"];
                (*octree)["pCount"] = (*octree)["pCount"].asUInt() + kv.second->size();
                for (size_t i=0; i<kv.first.length(); i++) {
                    octree = &(*octree)[kv.first.substr(i,1)];
                    if (*octree != Json::nullValue) {
                        (*octree)["pCount"] = (*octree)["pCount"].asUInt() + kv.second->size();
                    }
                }
            }
            writeIndexInfo();
            std::cout << "[" << currentTime() << "] add vertex data end." << std::endl;
        }

        void queryBoundingBox(const PCT::BoundingBox& box, PCT::PointCloudData& data) {
            Json::Value& octree = indexInfo["octree"];

            if (octree["pCount"].asUInt() <= 0) {
                return;
            }

            std::cout << "[" << currentTime() << "] query start, box{" << box.min[0] << "," << box.min[1] << ","  << box.min[2] << "-->" << box.max[0] << ","  << box.max[1] << ","  << box.max[2] << "}" << std::endl;

            data.reset();
            doQueryPolyhedron(box, data, octree);
            std::cout << "[" << currentTime() << "] query end." << std::endl;
        }

        void queryVisionCone(const PCT::VisionCone& visionCone, PCT::PointCloudData& data) {
            Json::Value& octree = indexInfo["octree"];

            if (octree["pCount"].asUInt() <= 0) {
                return;
            }

            std::cout << "[" << currentTime() << "] query start, perspective{fov:" << visionCone.fov << ", aspect:" << visionCone.aspect << ", near:"  << visionCone.near << ", far:" << visionCone.far << "}" << std::endl;

            data.reset();
            doQueryPolyhedron(visionCone, data, octree);
            std::cout << "[" << currentTime() << "] query end." << std::endl;
        }

        void queryData(PCT::PointCloudData& data, const Eigen::Matrix4d& viewMatrix, const float near, const float far, const float fov, const float aspect) {
            queryVisionCone(PCT::VisionCone(viewMatrix, near, far, fov, aspect), data);
        }

    private:
        void doQueryPolyhedron(const PCT::Polyhedron& polyhedron, PCT::PointCloudData& data, const Json::Value& octree) {
            if (octree["pCount"].asUInt() <= 0) {
                return;
            }
            int containType = polyhedron.checkContainment(octree["minX"].asDouble(), octree["minY"].asDouble(), octree["minZ"].asDouble(), octree["maxX"].asDouble(), octree["maxY"].asDouble(), octree["maxZ"].asDouble());
            if (containType == PCT::BoundingBox::RELATIONSHIP_NOT_INCLUDED) {
                return;
            }
            if (octree["depth"].asInt() == 0) {
                try {
                    std::cout << "[" << currentTime() << "] get PCD:" << octree["tag"].asString() << std::endl;
                    PCT::PLYData plyIn((cachePath/(octree["tag"].asString() + ".ply")).string());
                    plyIn.validate();
                    PCT::PointCloudData subData;
                    plyIn.getVertexData(subData);
                    if (data.getSpec().empty()) {
                        std::cout << "[" << currentTime() << "] get spec of PCD." << std::endl;
                        data.getSpec().insert(data.getSpec().end(), subData.getSpec().begin(), subData.getSpec().end());
                        data.notifySpecChanged();
                        data.initContainersBySpec();
                    }

                    std::cout << "[" << currentTime() << "] get data of PCD start." << std::endl;
                    if (containType == PCT::BoundingBox::RELATIONSHIP_FULL_INCLUDED) {
                        for (int i=0; i<data.getSpec().size(); i++) {
                            const std::pair<std::string,std::string>& pair = data.getSpec()[i];
                            const std::string& typeName = pair.first;
                            const std::string& typeStr = pair.second;
                            if (typeStr == "uchar" || typeStr == "uint8") {
                                pickDataToGeneralTypeVector<uint8_t>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "ushort" || typeStr == "uint16") {
                                pickDataToGeneralTypeVector<uint16_t>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "uint" || typeStr == "uint32") {
                                pickDataToGeneralTypeVector<uint32_t>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "char" || typeStr == "int8") {
                                pickDataToGeneralTypeVector<int8_t>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "short" || typeStr == "int16") {
                                pickDataToGeneralTypeVector<int16_t>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "int" || typeStr == "int32") {
                                pickDataToGeneralTypeVector<int32_t>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "float" || typeStr == "float32") {
                                pickDataToGeneralTypeVector<float>(data.getData()[i], subData.getData()[i]);
                            } else if (typeStr == "double" || typeStr == "float64") {
                                pickDataToGeneralTypeVector<double>(data.getData()[i], subData.getData()[i]);
                            }
                        }
                    } else if (containType == PCT::BoundingBox::RELATIONSHIP_PARTIALLY_INCLUDED) {
                        std::vector<size_t> indexes;
                        size_t count = subData.count();
                        for (size_t i=0; i<count; i++) {
                            if (polyhedron.ifContain(subData.x(i), subData.y(i), subData.z(i))) {
                                indexes.emplace_back(i);
                            }
                        }

                        if (!indexes.empty()) {
                            for (int i=0; i<data.getSpec().size(); i++) {
                                const std::pair<std::string,std::string>& pair = data.getSpec()[i];
                                const std::string& typeName = pair.first;
                                const std::string& typeStr = pair.second;
                                if (typeStr == "uchar" || typeStr == "uint8") {
                                    pickDataToGeneralTypeVector<uint8_t>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "ushort" || typeStr == "uint16") {
                                    pickDataToGeneralTypeVector<uint16_t>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "uint" || typeStr == "uint32") {
                                    pickDataToGeneralTypeVector<uint32_t>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "char" || typeStr == "int8") {
                                    pickDataToGeneralTypeVector<int8_t>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "short" || typeStr == "int16") {
                                    pickDataToGeneralTypeVector<int16_t>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "int" || typeStr == "int32") {
                                    pickDataToGeneralTypeVector<int32_t>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "float" || typeStr == "float32") {
                                    pickDataToGeneralTypeVector<float>(data.getData()[i], subData.getData()[i], indexes);
                                } else if (typeStr == "double" || typeStr == "float64") {
                                    pickDataToGeneralTypeVector<double>(data.getData()[i], subData.getData()[i], indexes);
                                }
                            }
                        }
                    }
                    std::cout << "[" << currentTime() << "] get data of PCD end." << std::endl;
                } catch (std::exception& exception) {
                    std::cout << "err:" << exception.what() << std::endl;
                }
            } else {
                for (int i=0; i<8; i++) {
                    doQueryPolyhedron(polyhedron, data, octree[std::to_string(i)]);
                }
            }
        }
    private:
        fs::path cachePath;
        fs::path pathIndexInfo;
        int depth{};
        PCT::BoundingBox boundingBox;
        Json::Value indexInfo;

        void initIndexInfo() {
            indexInfo = Json::Value();
            double maxLength = boundingBox.getMaxLength();
            indexInfo["octree"] = getOctree(depth,
                                            boundingBox.min[0],
                                            boundingBox.min[1],
                                            boundingBox.min[2],
                                            boundingBox.min[0] + maxLength,
                                            boundingBox.min[1] + maxLength,
                                            boundingBox.min[2] + maxLength,
                                            "");
            writeIndexInfo();
        }

        Json::Value getOctree(int octreeDepth, double minX, double minY, double minZ, double maxX, double maxY, double maxZ, const std::string& tag) {
            Json::Value value;
            value["minX"] = minX;
            value["minY"] = minY;
            value["minZ"] = minZ;
            value["maxX"] = maxX;
            value["maxY"] = maxY;
            value["maxZ"] = maxZ;
            value["depth"] = octreeDepth;
            value["pCount"] = 0;
            value["tag"] = tag;
            if (octreeDepth > 0) {
                int newDepth = octreeDepth -1;
                double centerX = (maxX+minX)/2;
                double centerY = (maxY+minY)/2;
                double centerZ = (maxZ+minZ)/2;

                //0 (min-center)
                value["0"] = getOctree(newDepth, minX, minY, minZ, centerX, centerY, centerZ, tag+"0");
                //1 (centerX_minY_minZ-maxX_centerY_centerZ)
                value["1"] = getOctree(newDepth, centerX, minY, minZ, maxX, centerY, centerZ, tag+"1");
                //2 (minX_centerY_minZ-centerX_maxY_centerZ)
                value["2"] = getOctree(newDepth, minX, centerY, minZ, centerX, maxY, centerZ, tag+"2");
                //3 (centerX_centerY_minZ-maxX_maxY_centerZ)
                value["3"] = getOctree(newDepth, centerX, centerY, minZ, maxX, maxY, centerZ, tag+"3");
                //4 (minX_minY_centerZ-centerX_centerY_maxZ)
                value["4"] = getOctree(newDepth, minX, minY, centerZ, centerX, centerY, maxZ, tag+"4");
                //5 (centerX_minY_centerZ-maxX_centerY_maxZ)
                value["5"] = getOctree(newDepth, centerX, minY, centerZ, maxX, centerY, maxZ, tag+"5");
                //6 (minX_centerY_centerZ-centerX_maxY_maxZ)
                value["6"] = getOctree(newDepth, minX, centerY, centerZ, centerX, maxY, maxZ, tag+"6");
                //7 (center-max)
                value["7"] = getOctree(newDepth, centerX, centerY, centerZ, maxX, maxY, maxZ, tag+"7");
            }
            return value;
        }

        void readIndexInfo() {
            if (fs::exists(pathIndexInfo)) {
                std::ifstream is(pathIndexInfo, std::ios::binary);
                if (!is.is_open()) {
                    throw std::runtime_error("[PCT]Open file failed:"+pathIndexInfo.string());
                }
                Json::Reader reader;
                if (reader.parse(is, indexInfo)) {
                    std::cout << "[PCT]Read octree success." << std::endl;
                    is.close();

                    this->depth = indexInfo["octree"]["depth"].asInt();
                } else {
                    is.close();
                    throw std::runtime_error("[PCT]Read octree failed.");
                }
            } else {
                throw std::runtime_error("[PCT]Index info is not exists!");
            }
        }

        void writeIndexInfo() {
            std::ofstream os(pathIndexInfo, std::ios::out);
            if (!os.is_open()) {
                throw std::runtime_error("[PCT]Open file failed:"+pathIndexInfo.string());
            }
            Json::StyledWriter sw;
            os << sw.write(indexInfo);
            os.close();
        }

        void writeSubPly(PCT::PointCloudData& data, const std::string& tag, const std::unique_ptr<std::vector<size_t>>& indexes) {
            const fs::path& path = cachePath / (tag + ".ply");
            PCT::PLYData plyOut;
            plyOut.addElement("vertex", indexes->size());
            PCT::Element& elemVertexOut = plyOut.getElement("vertex");
            for (int i=0; i<data.getSpec().size(); i++) {
                const std::pair<std::string,std::string>& pair = data.getSpec()[i];
                const std::string& typeName = pair.first;
                const std::string& typeStr = pair.second;
                if (typeStr == "uchar" || typeStr == "uint8") {
                    pickDataToElement<uint8_t>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "ushort" || typeStr == "uint16") {
                    pickDataToElement<uint16_t>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "uint" || typeStr == "uint32") {
                    pickDataToElement<uint32_t>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "char" || typeStr == "int8") {
                    pickDataToElement<int8_t>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "short" || typeStr == "int16") {
                    pickDataToElement<int16_t>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "int" || typeStr == "int32") {
                    pickDataToElement<int32_t>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "float" || typeStr == "float32") {
                    pickDataToElement<float>(typeName, data.getData()[i], indexes, elemVertexOut);
                } else if (typeStr == "double" || typeStr == "float64") {
                    pickDataToElement<double>(typeName, data.getData()[i], indexes, elemVertexOut);
                }
            }
            if (fs::exists(path)) {
                plyOut.append(path.string(), PCT::DataFormat::Binary);
            } else {
                plyOut.write(path.string(), PCT::DataFormat::Binary);
            }
        }

        template <typename T>
        void pickDataToElement(const std::string& typeName, std::shared_ptr<void> dataPtr, const std::unique_ptr<std::vector<size_t>>& indexes, PCT::Element& elemVertexOut) {
            auto propertyData = std::static_pointer_cast<std::vector<T>>(dataPtr);
            std::vector<T> vector;
            vector.reserve(indexes->size());
            for (const size_t index : *indexes) {
                vector.emplace_back((*propertyData)[index]);
            }
            elemVertexOut.addProperty(typeName, vector);
        }

        bool classify(const size_t pointIndex, const double& x, const double& y, const double& z, const Json::Value& octree, std::unordered_map<std::string,std::unique_ptr<std::vector<size_t>>>& octreeDataIndexMap) {
            if (x >= octree["minX"].asDouble() && x < octree["maxX"].asDouble()
                && y >= octree["minY"].asDouble() && y < octree["maxY"].asDouble()
                && z >= octree["minZ"].asDouble() && z < octree["maxZ"].asDouble()) {
                if (octree["depth"].asInt() <= 0) {
                    const std::string& tag = octree["tag"].asString();
                    auto got = octreeDataIndexMap.find(tag);
                    if (got == octreeDataIndexMap.end()) {
                        octreeDataIndexMap.emplace(tag, std::make_unique<std::vector<size_t>>());
                        octreeDataIndexMap[tag]->emplace_back(pointIndex);
                    } else {
                        got->second->emplace_back(pointIndex);
                    }
                    return true;
                } else {
                    for (int i=0; i<8; i++) {
                        if (classify(pointIndex, x, y, z, octree[std::to_string(i)], octreeDataIndexMap)) {
                            return true;
                        }
                    }
                }
            } else {
                return false;
            }
        }

        void fastClassify(const PCT::PointCloudData& data, const Json::Value& octree, std::unordered_map<std::string,std::unique_ptr<std::vector<size_t>>>& octreeDataIndexMap) {
            const double& originMinX = octree["minX"].asDouble();
            const double& originMinY = octree["minY"].asDouble();
            const double& originMinZ = octree["minZ"].asDouble();
            const double& originMaxX = octree["maxX"].asDouble();
            const double& originMaxY = octree["maxY"].asDouble();
            const double& originMaxZ = octree["maxZ"].asDouble();
            const int& originDp = octree["depth"].asInt();

            size_t count = data.count();
            for (size_t pointIndex = 0; pointIndex < count; pointIndex++) {
                const double& x = data.x(pointIndex);
                const double& y = data.y(pointIndex);
                const double& z = data.z(pointIndex);
                double minX = originMinX;
                double minY = originMinY;
                double minZ = originMinZ;
                double maxX = originMaxX;
                double maxY = originMaxY;
                double maxZ = originMaxZ;
                if (x >= minX && x < maxX && y >= minY && y < maxY && z >= minZ && z < maxZ) {
                    std::string tag;
                    for (int dp = originDp;dp > 0;dp--) {
                        double middleX = (maxX + minX)/2;
                        double middleY = (maxY + minY)/2;
                        double middleZ = (maxZ + minZ)/2;
                        bool halfX = x >= middleX;
                        bool halfY = y >= middleY;
                        bool halfZ = z >= middleZ;
                        std::string subTag;
                        if (!halfX && !halfY && !halfZ) {
                            tag += "0";
//                            minX = minX;
//                            minY = minY;
//                            minZ = minZ;
                            maxX = middleX;
                            maxY = middleY;
                            maxZ = middleZ;
                        } else if (halfX && !halfY && !halfZ) {
                            tag += "1";
                            minX = middleX;
//                            minY = minY;
//                            minZ = minZ;
//                            maxX = maxX;
                            maxY = middleY;
                            maxZ = middleZ;
                        } else if (!halfX && halfY && !halfZ) {
                            tag += "2";
//                            minX = minX;
                            minY = middleY;
//                            minZ = minZ;
                            maxX = middleX;
//                            maxY = maxY;
                            maxZ = middleZ;
                        } else if (halfX && halfY && !halfZ) {
                            tag += "3";
                            minX = middleX;
                            minY = middleY;
//                            minZ = minZ;
//                            maxX = maxX;
//                            maxY = maxY;
                            maxZ = middleZ;
                        } else if (!halfX && !halfY && halfZ) {
                            tag += "4";
//                            minX = minX;
//                            minY = minY;
                            minZ = middleZ;
                            maxX = middleX;
                            maxY = middleY;
//                            maxZ = maxZ;
                        } else if (halfX && !halfY && halfZ) {
                            tag += "5";
                            minX = middleX;
//                            minY = minY;
                            minZ = middleZ;
//                            maxX = maxX;
                            maxY = middleY;
//                            maxZ = maxZ;
                        } else if (!halfX && halfY && halfZ) {
                            tag += "6";
//                            minX = minX;
                            minY = middleY;
                            minZ = middleZ;
                            maxX = middleX;
//                            maxY = maxY;
//                            maxZ = maxZ;
                        } else if (halfX && halfY && halfZ) {
                            tag += "7";
                            minX = middleX;
                            minY = middleY;
                            minZ = middleZ;
//                            maxX = maxX;
//                            maxY = maxY;
//                            maxZ = maxZ;
                        }
                    }
                    if (tag.length() == 4) {
                        auto got = octreeDataIndexMap.find(tag);
                        if (got == octreeDataIndexMap.end()) {
                            octreeDataIndexMap.emplace(tag, std::make_unique<std::vector<size_t>>());
                            octreeDataIndexMap[tag]->emplace_back(pointIndex);
                        } else {
                            got->second->emplace_back(pointIndex);
                        }
                    } else {
                        throw std::runtime_error("unknown errors occur in fastClassify");
                    }
                }
            }

        }

        template<typename T>
        void pickDataToGeneralTypeVector(std::shared_ptr<void>& dest, const std::shared_ptr<void>& src) {
            std::shared_ptr<std::vector<T>> d = std::static_pointer_cast<std::vector<T>>(dest);
            std::shared_ptr<std::vector<T>> s = std::static_pointer_cast<std::vector<T>>(src);
            d->insert(d->end(), s->begin(), s->end());
        }

        template<typename T>
        void pickDataToGeneralTypeVector(std::shared_ptr<void>& dest, const std::shared_ptr<void>& src, const std::vector<size_t>& indexes) {
            std::shared_ptr<std::vector<T>> d = std::static_pointer_cast<std::vector<T>>(dest);
            std::shared_ptr<std::vector<T>> s = std::static_pointer_cast<std::vector<T>>(src);
            d->reserve(d->capacity() + indexes.size());
            for (const size_t& index : indexes) {
                d->emplace_back((*s)[index]);
            }
        }

        inline static std::string currentTime() {
            return currentFormatTime();
        }

        inline static time_t currentTimeMilliseconds() {
            return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        }

        inline static std::string currentFormatTime() {
            time_t t = time(nullptr);
            char buff[32]={NULL};
            tm ltm{0};
            localtime_s(&ltm, &t);
            strftime(buff, sizeof (buff), "%Y-%m-%d %H:%M:%S", &ltm);
            return buff;
        }
    };
}
#endif //PLYFASTLOADER_POINTCLOUNDSEARCHER_H
