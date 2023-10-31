//
// Created by Zachary on 2023/9/27.
//

#ifndef PLYFASTLOADER_PCT_ENTITY_H
#define PLYFASTLOADER_PCT_ENTITY_H
#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <array>
#include <eigen3/Eigen/Dense>

namespace PCT {
#define PI 3.141592654
    /**
     * An interface for multi-dimension data access.
     */
    class IMultiDimensionDataSet {
    public:
        virtual double dimensionValueAtIndex(size_t index, size_t dimension) const = 0;
        virtual size_t dataCount() const = 0;
    };

    class Polyhedron {
    public:
        virtual int checkContainment(const double& minX, const double& minY, const double& minZ, const double& maxX, const double& maxY, const double& maxZ) const = 0;
        virtual bool ifContain(const double& x, const double& y, const double& z) const = 0;
    };
    /**
     * min:inclusive, max:exclusive
     */
    class BoundingBox : public Polyhedron {

    public:
        const static int RELATIONSHIP_FULL_INCLUDED = 0;
        const static int RELATIONSHIP_PARTIALLY_INCLUDED = 1;
        const static int RELATIONSHIP_NOT_INCLUDED = 2;

        std::array<double, 3> min{};
        std::array<double, 3> max{};

        BoundingBox() = default;
        BoundingBox(double minX, double minY, double minZ, double maxX, double maxY, double maxZ) : min({minX, minY, minZ}), max({maxX, maxY, maxZ}) {}
        ~BoundingBox()= default;

        inline int checkContainment(const double &minX, const double &minY, const double &minZ, const double &maxX,
                             const double &maxY, const double &maxZ) const override {
            if (minX > this->min[0] && minY > this->min[1] && minZ > this->min[2] &&
                maxX < this->max[0] && maxY < this->max[1] && maxZ < this->max[2]) {
                return RELATIONSHIP_FULL_INCLUDED;
            } else if (this->min[0] > maxX || this->min[1] > maxY || this->min[2] > maxZ ||
                       this->max[0] < minX || this->max[1] < minY || this->max[2] < minZ) {
                return RELATIONSHIP_NOT_INCLUDED;
            } else {
                return RELATIONSHIP_PARTIALLY_INCLUDED;
            }
        }

        inline bool ifContain(const double &x, const double &y, const double &z) const override {
            return x >= this->min[0] && y >= this->min[1] && z >=this->min[2] &&
                   x < this->max[0] && y < this->max[1] && z < this->max[2];
        }

        inline double getMaxLength() const {
            double a = max[0]-min[0], b = max[1]-min[1],c = max[2]-min[2];
            return std::max(std::max(a,b),c);
        }
    };

    class VisionCone : public Polyhedron {

    public:
        const static int RELATIONSHIP_FULL_INCLUDED = 0;
        const static int RELATIONSHIP_PARTIALLY_INCLUDED = 1;
        const static int RELATIONSHIP_NOT_INCLUDED = 2;

        const float near, far, fov, aspect;
        const Eigen::Matrix4d vMatrix;

    private:
        Eigen::Matrix4d pMatrix;
        Eigen::Matrix4d pvMatrix;

        std::array<double,3> min{-1, -1, -1}, max{1, 1, 1};//标准设备空间范围

    public:
        VisionCone(Eigen::Matrix4d viewMatrix, float near, float far, float fov, float aspect) : vMatrix(std::move(viewMatrix)), near(near), far(far), fov(fov), aspect(aspect) {
            perspectiveM(pMatrix, fov, aspect, near, far);
            pvMatrix = pMatrix*vMatrix;
        }
        ~VisionCone()= default;

        int checkContainment(const double &minX, const double &minY, const double &minZ, const double &maxX,
                             const double &maxY, const double &maxZ) const override {
//            std::cout << "checkContainment,Box-{"<<minX<<","<<minY<<","<<minZ<<"}{"<<maxX<<","<<maxY<<","<<maxZ<<"}"<<std::endl;
            std::array<Eigen::Vector4d, 8> vertexes;
            //待测包围盒各顶点变换到观察空间
            vertexes[0] = {minX,minY,minZ,1};
            vertexes[1] = {maxX,minY,minZ,1};
            vertexes[2] = {minX,maxY,minZ,1};
            vertexes[3] = {maxX,maxY,minZ,1};
            vertexes[4] = {minX,minY,maxZ,1};
            vertexes[5] = {maxX,minY,maxZ,1};
            vertexes[6] = {minX,maxY,maxZ,1};
            vertexes[7] = {maxX,maxY,maxZ,1};
            for (Eigen::Vector4d& vertex : vertexes) {
                //变换到观察空间
                vertex = vMatrix*vertex;
            }
            //获取在观察空间的新包围盒
            double tMinX = vertexes[0](0), tMinY = vertexes[0](1), tMinZ = vertexes[0](2), tMaxX = vertexes[0](0), tMaxY = vertexes[0](1), tMaxZ = vertexes[0](2);
            for (Eigen::Vector4d& vertex : vertexes) {
                if (vertex(0) < tMinX) {
                    tMinX = vertex(0);
                }
                if (vertex(0) > tMaxX) {
                    tMaxX = vertex(0);
                }
                if (vertex(1) < tMinY) {
                    tMinY = vertex(1);
                }
                if (vertex(1) > tMaxY) {
                    tMaxY = vertex(1);
                }
                if (vertex(2) < tMinZ) {
                    tMinZ = vertex(2);
                }
                if (vertex(2) > tMaxZ) {
                    tMaxZ = vertex(2);
                }
            }
//            std::cout << "checkContainment,tBox-{"<<tMinX<<","<<tMinY<<","<<tMinZ<<"}{"<<tMaxX<<","<<tMaxY<<","<<tMaxZ<<"}"<<std::endl;

            if (tMaxZ < -far || tMinZ > -near) {
                return RELATIONSHIP_NOT_INCLUDED;
            }

            double tanYZ = std::tan(fov * (PI / 360.0));
            double tanXZ = tanYZ*aspect;
            double lz = std::abs(tMinZ);
            if ((tMinY > 0 && tMinY/lz > tanYZ) || (tMaxY < 0 && std::abs(tMaxY)/lz > tanYZ) || (tMinX > 0 && tMinX/lz > tanXZ) || (tMaxX < 0 && std::abs(tMaxX)/lz > tanXZ)) {
                return RELATIONSHIP_NOT_INCLUDED;
            }

            lz = std::abs(tMaxZ);
            if (tMinZ > -far && tMaxZ < -near && std::abs(tMaxY)/lz <= tanYZ && std::abs(tMinY)/lz <= tanYZ && std::abs(tMaxX)/lz <= tanXZ && std::abs(tMinX)/lz <= tanXZ) {
                return RELATIONSHIP_FULL_INCLUDED;
            }

            return RELATIONSHIP_PARTIALLY_INCLUDED;
        }

        inline bool ifContain(const double &x, const double &y, const double &z) const override {
            Eigen::Vector4d p(x, y, z, 1);
            //变换到投影空间
            p = pvMatrix*p;
            //归一化到标准设备空间
            p /= p(3);

            return p(0) >= this->min[0] && p(1) >= this->min[1] && p(2) >=this->min[2] &&
                   p(0) <= this->max[0] && p(1) <= this->max[1] && p(2) <= this->max[2];
        }
    private:
        static void perspectiveM(Eigen::Matrix4d& m, float fovy, float aspect, float zNear, float zFar) {
            float f = 1.0f / (float) std::tan(fovy * (PI / 360.0));
            float rangeReciprocal = 1.0f / (zNear - zFar);
            m << f / aspect, 0.0f, 0.0f, 0.0f,
                    0.0f, f, 0.0f, 0.0f,
                    0.0f, 0.0f, (zFar + zNear) * rangeReciprocal, 2.0f * zFar * zNear * rangeReciprocal,
                    0.0f, 0.0f, -1.0f, 0.0f;
        }
    };

    /**
     * vertex element data in ply
     */
    class PointCloudData : public IMultiDimensionDataSet {
    public:
        PointCloudData() = default;
        PointCloudData(const std::vector<std::pair<std::string,std::string>>& spec, const std::vector<std::shared_ptr<void>>& data) {
            this->spec = spec;
            this->data = data;
            notifySpecChanged();
            notifyDataContainersChanged();
        }
        ~PointCloudData() = default;

        inline std::vector<std::pair<std::string, std::string>> &getSpec() {
            return spec;
        }

        inline std::vector<std::shared_ptr<void>> &getData() {
            return data;
        }

        inline const std::vector<std::shared_ptr<void>> &getData() const {
            return data;
        }

        inline void reset() {
            spec.clear();
            data.clear();
        }

        void initContainersBySpec() {
            initContainersBySpec(data);
        }

        void notifySpecChanged() {
            for (int i=0; i<spec.size(); i++) {
                const std::pair<std::string,std::string>& pair = spec[i];
                if (pair.first == "x") {
                    dimensionIndexes[0] = i;
                    dimensionTypes[0] = getValueType(pair.second);
                } else if (pair.first == "y") {
                    dimensionIndexes[1] = i;
                    dimensionTypes[1] = getValueType(pair.second);
                } else if (pair.first == "z") {
                    dimensionIndexes[2] = i;
                    dimensionTypes[2] = getValueType(pair.second);
                } else if (pair.first == "red") {
                    rIndex = i;
                    valueTypeR = getValueType(pair.second);
                } else if (pair.first == "green") {
                    gIndex = i;
                    valueTypeG = getValueType(pair.second);
                } else if (pair.first == "blue") {
                    bIndex = i;
                    valueTypeB = getValueType(pair.second);
                } else if (pair.first == "intensity") {
                    intensityIndex = i;
                    valueTypeIntensity = getValueType(pair.second);
                }
            }
        }

        void notifyDataContainersChanged() {
            for (int dimension=0; dimension<dimensionTypes.size(); dimension++) {
                if (dimensionTypes[dimension] == VALUE_TYPE_FLOAT) {
                    floatDimensionDataSharedPtr[dimension] = std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[dimension]]);
                    floatDimensionDataPtr[dimension] = floatDimensionDataSharedPtr[dimension].get();
                } else if (dimensionTypes[dimension] == VALUE_TYPE_DOUBLE) {
                    doubleDimensionDataSharedPtr[dimension] = std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[dimension]]);
                    doubleDimensionDataPtr[dimension] = doubleDimensionDataSharedPtr[dimension].get();
                } else {
                    throw std::runtime_error("Invalid dimension value type:" + std::to_string(dimensionTypes[dimension]));
                }
            }
        }

        inline size_t count() const {
            if (dimensionTypes[0] == VALUE_TYPE_FLOAT) {
                return std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[0]])->size();
            } else if (dimensionTypes[0] == VALUE_TYPE_DOUBLE) {
                return std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[0]])->size();
            } else if (spec.empty() || data.empty()) {
                return 0;
            } else {
                throw std::runtime_error("Invalid x value type in parsing points count.");
            }
        }

        double dimensionValueAtIndex(size_t index, size_t dimension) const override {

            if (dimensionTypes[dimension] == VALUE_TYPE_FLOAT) {
                return (*floatDimensionDataPtr[dimension])[index];
//                return (*std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[dimension]]))[index];
            } else if (dimensionTypes[dimension] == VALUE_TYPE_DOUBLE) {
                return (*doubleDimensionDataPtr[dimension])[index];
//                return (*std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[dimension]]))[index];
            } else {
                throw std::runtime_error("Invalid dimension value type:" + std::to_string(dimensionTypes[dimension]));
            }
        }

        size_t dataCount() const override {
            return count();
        }

        inline double x(size_t index) const {
            if (dimensionTypes[0] == VALUE_TYPE_FLOAT) {
                return std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[0]])->at(index);
            } else if (dimensionTypes[0] == VALUE_TYPE_DOUBLE) {
                return std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[0]])->at(index);
            } else {
                throw std::runtime_error("Invalid x value type:" + std::to_string(dimensionTypes[0]));
            }
        }

        inline double y(size_t index) const {
            if (dimensionTypes[1] == VALUE_TYPE_FLOAT) {
                return std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[1]])->at(index);
            } else if (dimensionTypes[1] == VALUE_TYPE_DOUBLE) {
                return std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[1]])->at(index);
            } else {
                throw std::runtime_error("Invalid y value type:" + std::to_string(dimensionTypes[1]));
            }
        }

        inline double z(size_t index) const {
            if (dimensionTypes[2] == VALUE_TYPE_FLOAT) {
                return std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[2]])->at(index);
            } else if (dimensionTypes[2] == VALUE_TYPE_DOUBLE) {
                return std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[2]])->at(index);
            } else {
                throw std::runtime_error("Invalid z value type:" + std::to_string(dimensionTypes[2]));
            }
        }

        inline uint8_t r(size_t index) const {
            if (valueTypeR == VALUE_TYPE_UINT8) {
                return std::static_pointer_cast<std::vector<uint8_t>>(data[rIndex])->at(index);
            } else {
                throw std::runtime_error("Invalid r value type:" + std::to_string(valueTypeR));
            }
        }

        inline uint8_t g(size_t index) const {
            if (valueTypeG == VALUE_TYPE_UINT8) {
                return std::static_pointer_cast<std::vector<uint8_t>>(data[gIndex])->at(index);
            } else {
                throw std::runtime_error("Invalid g value type:" + std::to_string(valueTypeG));
            }
        }

        inline uint8_t b(size_t index) const {
            if (valueTypeB == VALUE_TYPE_UINT8) {
                return std::static_pointer_cast<std::vector<uint8_t>>(data[bIndex])->at(index);
            } else {
                throw std::runtime_error("Invalid b value type:" + std::to_string(valueTypeB));
            }
        }

        inline float intensity(size_t index) const {
            if (valueTypeIntensity == VALUE_TYPE_FLOAT) {
                return std::static_pointer_cast<std::vector<float>>(data[intensityIndex])->at(index);
            }  else {
                throw std::runtime_error("Invalid z value type:" + std::to_string(valueTypeIntensity));
            }
        }

        inline void setX(size_t index, const double value) {
            if (dimensionTypes[0] == VALUE_TYPE_FLOAT) {
                auto& v= std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[0]])->at(index);
                v = static_cast<float> (value);
            } else if (dimensionTypes[0] == VALUE_TYPE_DOUBLE) {
                auto& v= std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[0]])->at(index);
                v = value;
            } else {
                throw std::runtime_error("Invalid x value type.");
            }
        }

        inline void setY(size_t index, const double value) {
            if (dimensionTypes[1] == VALUE_TYPE_FLOAT) {
                auto& v= std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[1]])->at(index);
                v = static_cast<float> (value);
            } else if (dimensionTypes[1] == VALUE_TYPE_DOUBLE) {
                auto& v= std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[1]])->at(index);
                v = value;
            } else {
                throw std::runtime_error("Invalid y value type.");
            }
        }

        inline void setZ(size_t index, const double value) {
            if (dimensionTypes[2] == VALUE_TYPE_FLOAT) {
                auto& v= std::static_pointer_cast<std::vector<float>>(data[dimensionIndexes[2]])->at(index);
                v = static_cast<float> (value);
            } else if (dimensionTypes[2] == VALUE_TYPE_DOUBLE) {
                auto& v= std::static_pointer_cast<std::vector<double>>(data[dimensionIndexes[2]])->at(index);
                v = value;
            } else {
                throw std::runtime_error("Invalid z value type.");
            }
        }

        void transform(const Eigen::Matrix4d& mt) {
            size_t count = this->count();
            double p[3];
            for (size_t i=0; i<count; i++) {
                p[0] = x(i);
                p[1] = y(i);
                p[2] = z(i);
                setX(i, mt(0, 0) * p[0] + mt(0, 1) * p[1] + mt(0, 2) * p[2] + mt(0, 3));
                setY(i, mt(1, 0) * p[0] + mt(1, 1) * p[1] + mt(1, 2) * p[2] + mt(1, 3));
                setZ(i, mt(2, 0) * p[0] + mt(2, 1) * p[1] + mt(2, 2) * p[2] + mt(2, 3));
            }
        }

        void crop(const BoundingBox& box) {
            std::vector<size_t> indexes;
            size_t count = this->count();
            for (size_t i=0; i<count; i++) {
                if (box.ifContain(x(i), y(i), z(i))) {
                    indexes.emplace_back(i);
                }
            }

            std::vector<std::shared_ptr<void>> newData;
            initContainersBySpec(newData);
            for (int i=0; i<spec.size(); i++) {
                const std::string& typeStr = spec[i].second;
                copy(data[i], newData[i], indexes, typeStr);
            }
            data.clear();
            data.insert(data.end(), newData.begin(), newData.end());
        }

        void subCloud(PointCloudData& target, std::vector<size_t>& indexes) {
            if (!target.getSpec().empty() || !target.getData().empty()) {
                throw std::runtime_error("The target sub cloud data is not empty.");
            }
            target.getSpec().insert(target.getSpec().end(), spec.begin(), spec.end());
            target.notifySpecChanged();
            target.initContainersBySpec();

            if (!indexes.empty()) {
                for (int i=0; i<spec.size(); i++) {
                    const std::string& typeStr = spec[i].second;
                    if (typeStr == "uchar" || typeStr == "uint8") {
                        pickDataToGeneralTypeVector<uint8_t>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "ushort" || typeStr == "uint16") {
                        pickDataToGeneralTypeVector<uint16_t>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "uint" || typeStr == "uint32") {
                        pickDataToGeneralTypeVector<uint32_t>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "char" || typeStr == "int8") {
                        pickDataToGeneralTypeVector<int8_t>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "short" || typeStr == "int16") {
                        pickDataToGeneralTypeVector<int16_t>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "int" || typeStr == "int32") {
                        pickDataToGeneralTypeVector<int32_t>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "float" || typeStr == "float32") {
                        pickDataToGeneralTypeVector<float>(data[i], target.getData()[i], indexes);
                    } else if (typeStr == "double" || typeStr == "float64") {
                        pickDataToGeneralTypeVector<double>(data[i], target.getData()[i], indexes);
                    }
                }
            }
        }

        void toNDC(const Eigen::Matrix4d& viewMatrix, float near, float far, float fov, float aspect) {
            Eigen::Matrix4d pMatrix;
            perspectiveM(pMatrix, fov, aspect, near, far);
            Eigen::Matrix4d pvMatrix = pMatrix*viewMatrix;
            const size_t count = dataCount();
            for (int i=0; i<count; i++) {
                Eigen::Vector4d p(x(i), y(i), z(i), 1);
                //变换到投影空间
                p = pvMatrix*p;
                //归一化到标准设备空间
                p /= p(3);

                setX(i, p(0));
                setY(i, p(1));
                setZ(i, p(2));
            }
        }

        /**
         * operator +=
         * @param other
         * @return this reference
         */
        PointCloudData& operator+=(const PointCloudData& other) {
            if (other.spec.size() != this->spec.size()) {
                throw std::runtime_error("Invalid spec size");
            }
            for (int i=0; i<other.spec.size(); i++) {
                if (other.spec[i] != this->spec[i]) {
                    throw std::runtime_error("Invalid spec, need{" + this->spec[i].first + ":" + this->spec[i].second + "} but {" + other.spec[i].first + ":" + other.spec[i].second + "} provided.");
                }
            }
            for (int i=0; i<spec.size(); i++) {
                const std::string& typeStr = spec[i].second;
                if (typeStr == "uchar" || typeStr == "uint8") {
                    pickDataToGeneralTypeVector<uint8_t>(other.getData()[i], this->data[i]);
                } else if (typeStr == "ushort" || typeStr == "uint16") {
                    pickDataToGeneralTypeVector<uint16_t>(other.getData()[i], this->data[i]);
                } else if (typeStr == "uint" || typeStr == "uint32") {
                    pickDataToGeneralTypeVector<uint32_t>(other.getData()[i], this->data[i]);
                } else if (typeStr == "char" || typeStr == "int8") {
                    pickDataToGeneralTypeVector<int8_t>(other.getData()[i], this->data[i]);
                } else if (typeStr == "short" || typeStr == "int16") {
                    pickDataToGeneralTypeVector<int16_t>(other.getData()[i], this->data[i]);
                } else if (typeStr == "int" || typeStr == "int32") {
                    pickDataToGeneralTypeVector<int32_t>(other.getData()[i], this->data[i]);
                } else if (typeStr == "float" || typeStr == "float32") {
                    pickDataToGeneralTypeVector<float>(other.getData()[i], this->data[i]);
                } else if (typeStr == "double" || typeStr == "float64") {
                    pickDataToGeneralTypeVector<double>(other.getData()[i], this->data[i]);
                }
            }
            return *this;
        }

    private:
        std::vector<std::pair<std::string,std::string>> spec;
        std::vector<std::shared_ptr<void>> data;

        std::array<std::shared_ptr<std::vector<float>>, 3> floatDimensionDataSharedPtr;
        std::array<std::shared_ptr<std::vector<double>>, 3> doubleDimensionDataSharedPtr;
        std::array<std::vector<float>*, 3> floatDimensionDataPtr;
        std::array<std::vector<double>*,3> doubleDimensionDataPtr;

        std::array<size_t,3> dimensionIndexes{};
        std::array<int,3> dimensionTypes{};

        const static int VALUE_TYPE_INVALID = -1, VALUE_TYPE_UINT8 = 0, VALUE_TYPE_UNIT16 = 1, VALUE_TYPE_UINT32 = 2, VALUE_TYPE_INT8 = 3, VALUE_TYPE_INT16 = 4, VALUE_TYPE_INT32 = 5, VALUE_TYPE_FLOAT = 6, VALUE_TYPE_DOUBLE = 7;
        int rIndex = -1, gIndex = -1, bIndex = -1, intensityIndex = -1;
        int valueTypeR = VALUE_TYPE_INVALID;
        int valueTypeG = VALUE_TYPE_INVALID;
        int valueTypeB = VALUE_TYPE_INVALID;
        int valueTypeIntensity = VALUE_TYPE_INVALID;

        static int getValueType(const std::string& typeStr) {
            if (typeStr == "float" || typeStr == "float32") {
                return VALUE_TYPE_FLOAT;
            } else if (typeStr == "double" || typeStr == "float64") {
                return VALUE_TYPE_DOUBLE;
            } else if (typeStr == "uchar" || typeStr == "uint8") {
                return VALUE_TYPE_UINT8;
            } else if (typeStr == "ushort" || typeStr == "uint16") {
                return VALUE_TYPE_UNIT16;
            } else if (typeStr == "uint" || typeStr == "uint32") {
                return VALUE_TYPE_UINT32;
            } else if (typeStr == "char" || typeStr == "int8") {
                return VALUE_TYPE_INT8;
            } else if (typeStr == "short" || typeStr == "int16") {
                return VALUE_TYPE_INT16;
            } else if (typeStr == "int" || typeStr == "int32") {
                return VALUE_TYPE_INT32;
            } else {
                return VALUE_TYPE_INVALID;
            }
        }

        inline static void copy(const std::shared_ptr<void>& src, std::shared_ptr<void>& dest, const std::string& typeStr) {
            if (typeStr == "uchar" || typeStr == "uint8") {
                copy<uint8_t>(src, dest);
            } else if (typeStr == "ushort" || typeStr == "uint16") {
                copy<uint16_t>(src, dest);
            } else if (typeStr == "uint" || typeStr == "uint32") {
                copy<uint32_t>(src, dest);
            } else if (typeStr == "char" || typeStr == "int8") {
                copy<int8_t>(src, dest);
            } else if (typeStr == "short" || typeStr == "int16") {
                copy<int16_t>(src, dest);
            } else if (typeStr == "int" || typeStr == "int32") {
                copy<int32_t>(src, dest);
            } else if (typeStr == "float" || typeStr == "float32") {
                copy<float>(src, dest);
            } else if (typeStr == "double" || typeStr == "float64") {
                copy<double>(src, dest);
            }
        }

        inline static void copy(const std::shared_ptr<void>& src, std::shared_ptr<void>& dest, const std::vector<size_t>& indexes, const std::string& typeStr) {
            if (typeStr == "uchar" || typeStr == "uint8") {
                copy<uint8_t>(src, dest, indexes);
            } else if (typeStr == "ushort" || typeStr == "uint16") {
                copy<uint16_t>(src, dest, indexes);
            } else if (typeStr == "uint" || typeStr == "uint32") {
                copy<uint32_t>(src, dest, indexes);
            } else if (typeStr == "char" || typeStr == "int8") {
                copy<int8_t>(src, dest, indexes);
            } else if (typeStr == "short" || typeStr == "int16") {
                copy<int16_t>(src, dest, indexes);
            } else if (typeStr == "int" || typeStr == "int32") {
                copy<int32_t>(src, dest, indexes);
            } else if (typeStr == "float" || typeStr == "float32") {
                copy<float>(src, dest, indexes);
            } else if (typeStr == "double" || typeStr == "float64") {
                copy<double>(src, dest, indexes);
            }
        }

        template<typename T>
        inline static void copy(const std::shared_ptr<void>& src, std::shared_ptr<void>& dest) {
            std::shared_ptr<std::vector<T>> d = std::static_pointer_cast<std::vector<T>>(dest);
            std::shared_ptr<std::vector<T>> s = std::static_pointer_cast<std::vector<T>>(src);
            d->insert(d->end(), s->begin(), s->end());
        }

        template<typename T>
        inline static void copy(const std::shared_ptr<void>& src, std::shared_ptr<void>& dest, const std::vector<size_t>& indexes) {
            std::shared_ptr<std::vector<T>> d = std::static_pointer_cast<std::vector<T>>(dest);
            std::shared_ptr<std::vector<T>> s = std::static_pointer_cast<std::vector<T>>(src);
            d->reserve(d->capacity() + indexes.size());
            for (const size_t& index : indexes) {
                d->emplace_back((*s)[index]);
            }
        }

        void initContainersBySpec(std::vector<std::shared_ptr<void>>& d) {
            d.clear();
            d.reserve(spec.size());
            for (int i=0; i<spec.size(); i++) {
                const std::string& typeStr = spec[i].second;
                if (typeStr == "uchar" || typeStr == "uint8") {
                    d.emplace_back(std::make_shared<std::vector<uint8_t>>());
                } else if (typeStr == "ushort" || typeStr == "uint16") {
                    d.emplace_back(std::make_shared<std::vector<uint16_t>>());
                } else if (typeStr == "uint" || typeStr == "uint32") {
                    d.emplace_back(std::make_shared<std::vector<uint32_t>>());
                } else if (typeStr == "char" || typeStr == "int8") {
                    d.emplace_back(std::make_shared<std::vector<int8_t>>());
                } else if (typeStr == "short" || typeStr == "int16") {
                    d.emplace_back(std::make_shared<std::vector<int16_t>>());
                } else if (typeStr == "int" || typeStr == "int32") {
                    d.emplace_back(std::make_shared<std::vector<int32_t>>());
                } else if (typeStr == "float" || typeStr == "float32") {
                    d.emplace_back(std::make_shared<std::vector<float>>());
                } else if (typeStr == "double" || typeStr == "float64") {
                    d.emplace_back(std::make_shared<std::vector<double>>());
                }
            }
            notifyDataContainersChanged();
        }

        template<typename T>
        inline void pickDataToGeneralTypeVector(const std::shared_ptr<void>& src, std::shared_ptr<void>& dest) {
            std::shared_ptr<std::vector<T>> s = std::static_pointer_cast<std::vector<T>>(src);
            std::shared_ptr<std::vector<T>> d = std::static_pointer_cast<std::vector<T>>(dest);
            d->insert(d->end(), s->begin(), s->end());
        }

        template<typename T>
        inline void pickDataToGeneralTypeVector(const std::shared_ptr<void>& src, std::shared_ptr<void>& dest, const std::vector<size_t>& indexes) {
            std::shared_ptr<std::vector<T>> s = std::static_pointer_cast<std::vector<T>>(src);
            std::shared_ptr<std::vector<T>> d = std::static_pointer_cast<std::vector<T>>(dest);
            d->reserve(d->capacity() + indexes.size());
            for (const size_t& index : indexes) {
                d->emplace_back((*s)[index]);
            }
        }

        static void perspectiveM(Eigen::Matrix4d& m, float fovy, float aspect, float zNear, float zFar) {
            float f = 1.0f / (float) std::tan(fovy * (PI / 360.0));
            float rangeReciprocal = 1.0f / (zNear - zFar);
            m << f / aspect, 0.0f, 0.0f, 0.0f,
                    0.0f, f, 0.0f, 0.0f,
                    0.0f, 0.0f, (zFar + zNear) * rangeReciprocal, 2.0f * zFar * zNear * rangeReciprocal,
                    0.0f, 0.0f, -1.0f, 0.0f;
        }

    public:
        std::string toString() const {
            std::string s = "PointCloudData->{count: " + std::to_string(count()) + "; spec: ";
            if (spec.size() > 0) {
                for (size_t i=0; i<spec.size(); i++) {
                    if (i > 0) {
                        s += ",";
                    }
                    s += spec[i].first + "[" + spec[i].second + "]";
                }
            }
            s += "}";
            return s;
        }
    };

    class TestData : public IMultiDimensionDataSet {
    public:
        TestData() {
            for (int i=0; i<2097152; i++) {
                float x = i;
                float y = i+1;
                float z = i+2;
                points.emplace_back(std::array<float,3>({x,y,z}));
            }
        }
        double dimensionValueAtIndex(size_t index, size_t dimension) const override {
            return points[index][dimension];
        }

        size_t dataCount() const override {
            return points.size();
        }
    private:
        std::vector<std::array<float, 3>> points;
    };
}
#endif //PLYFASTLOADER_PCT_ENTITY_H
