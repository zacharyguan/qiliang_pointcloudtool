//
// Created by Zachary on 2022/9/19.
//

#ifndef PLYFASTLOADER_TEST_H
#define PLYFASTLOADER_TEST_H
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "pointcloundsearcher.h"
#include "kdtree.h"

//#define CREATE
#define SEARCH
//#define CROP
//#define TRANSFORM
class Test {
public:
    Test() {
        using namespace PCT;
        const std::string& testFilePath = "E:\\local\\200.ply";//fulicaise
#if 0
        double minX=0,minY=0,minZ=0,maxX=0,maxY=0,maxZ=0;
        bool inited = false;
        size_t splitSize = 1024*1024*100;
        size_t blockCount = PLYData::prepareSplit(testFilePath, splitSize, false);
        std::cout << "blockCount:" << blockCount <<std::endl;
        if (blockCount > 0) {
            std::cout << "start calculate bounding box." << std::endl;
            for (int blockIndex=0; blockIndex<blockCount; blockIndex++) {
                PLYData plyIn(testFilePath, splitSize, blockIndex, false);
                try {
                    plyIn.validate();
                    PointCloudData pcd;
                    plyIn.getVertexData(pcd);
                    size_t count = pcd.count();
                    for (int i=0; i < count; i++) {
                        if (!inited) {
                            inited = true;
                            minX = maxX = pcd.x(0);
                            minY = maxY = pcd.y(0);
                            minZ = maxZ = pcd.z(0);
                        }
                        double v = pcd.x(i);
                        if (v < minX) {
                            minX = v;
                        }
                        if (v > maxX) {
                            maxX = v;
                        }
                        v = pcd.y(i);
                        if (v < minY) {
                            minY = v;
                        }
                        if (v > maxY) {
                            maxY = v;
                        }
                        v = pcd.z(i);
                        if (v < minZ) {
                            minZ = v;
                        }
                        if (v > maxZ) {
                            maxZ = v;
                        }
                    }
                } catch (std::exception& exception) {
                    std::cout << "err:" << exception.what() << std::endl;
                }
            }

        }
        std::cout << "BoundingBox-{"<<minX<<","<<minY<<","<<minZ<<"}{"<<maxX<<","<<maxY<<","<<maxZ<<"}"<<std::endl;
#endif
#ifdef CREATE
        if (blockCount > 0) {
            PCT::PointCloudSearcher pointCloudSearcher(R"(E:\local\pcs)", 4, BoundingBox(minX, minY, minZ, maxX, maxY, maxZ), true);

            for (int blockIndex=0; blockIndex<blockCount; blockIndex++) {
                PLYData plyIn(testFilePath, splitSize, blockIndex, false);
                try {
                    plyIn.validate();

                    PointCloudData tempData;
                    plyIn.getVertexData(tempData);

                    pointCloudSearcher.addVertexData(tempData);
                } catch (std::exception& exception) {
                    std::cout << "err:" << exception.what() << std::endl;
                }
            }
        }
#endif
#ifdef SEARCH
        PointCloudSearcher searcher(R"(E:\local\pcs\)");

        PointCloudData vertexData;

        Eigen::Matrix4d viewMatrix;
        viewMatrix.setIdentity();
        std::array<float, 16> temp{0};
        setLookAtM(temp, 0, 0, 0, 1, 0, 0, 0, 1, 0);
        viewMatrix << temp[0], temp[4], temp[8], temp[12],
                temp[1], temp[5], temp[9], temp[13],
                temp[2], temp[6], temp[10], temp[14],
                temp[3], temp[7], temp[11], temp[15];
        searcher.queryData(vertexData, viewMatrix, 0.001, 4, 20, 1);

//        searcher.queryBoundingBox(BoundingBox(minX, minY, minZ, maxX, maxY, maxZ), vertexData);

        for (int i=0; i<vertexData.getSpec().size(); i++) {
            const std::pair<std::string,std::string>& pair = vertexData.getSpec()[i];
            const std::string& typeName = pair.first;
            const std::string& typeStr = pair.second;
            if (typeStr == "uchar" || typeStr == "uint8") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<uint8_t>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "ushort" || typeStr == "uint16") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<uint16_t>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "uint" || typeStr == "uint32") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<uint32_t>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "char" || typeStr == "int8") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<int8_t>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "short" || typeStr == "int16") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<int16_t>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "int" || typeStr == "int32") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<int32_t>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "float" || typeStr == "float32") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<float>>(vertexData.getData()[i])->size() << std::endl;
            } else if (typeStr == "double" || typeStr == "float64") {
                std::cout << "Property " << typeName << ", count:" << std::static_pointer_cast<std::vector<double>>(vertexData.getData()[i])->size() << std::endl;
            }
        }
        PLYData::writeToPly(R"(E:\local\pcs\queryResult.ply)", vertexData);
#endif
#ifdef CROP
        PLYData plyIn(testFilePath);
        try {
            plyIn.validate();

            PointCloudData tempData;
            plyIn.getVertexData(tempData);
            tempData.crop(BoundingBox((minX+maxX)/2, minY, minZ, maxX, maxY, maxZ));
            PLYData::writeToPly(R"(E:\local\pcs\cropResult.ply)", tempData);
        } catch (std::exception& exception) {
            std::cout << "err:" << exception.what() << std::endl;
        }
#endif
#ifdef TRANSFORM
        {
            PLYData plyIn(testFilePath);
            try {
                plyIn.validate();

                PointCloudData tempData;
                plyIn.getVertexData(tempData);
                std::cout << "r:" << std::to_string(tempData.r(0)) << std::endl;
//                Eigen::Matrix4d mt;
//                mt << 0.996196925640,-0.048618465662,0.072304449975,-2.000000000000,
//                        0.052208472043,0.997452020645,-0.048618461937,1.000000000000,
//                        -0.069756470621,0.052208472043,0.996196925640,1.000000000000,
//                        0.000000000000,0.000000000000,0.000000000000,1.000000000000;
//                tempData.transform(mt);
//                PLYData::writeToPly(R"(E:\local\pcs\transformResult.ply)", tempData);
            } catch (std::exception& exception) {
                std::cout << "err:" << exception.what() << std::endl;
            }
        };
#endif
#if 1
        {
            try {
//                PLYData plyIn(testFilePath);
//                plyIn.validate();
//                PointCloudData tempData;
//                plyIn.getVertexData(tempData);
//                const size_t count = tempData.count();
//                std::cout << "read pcd size:" << std::to_string(count) << std::endl;
//                PointCloudData subData;
//                std::vector<size_t> indexes{0,1,2,3,4,5,6,7,8,9,10};
//                tempData.subCloud(subData, indexes);
//                std::cout << tempData.toString() <<std::endl;
//                std::cout << subData.toString() << std::endl;
////                    points.reserve(count);
////                    std::cout << "after reserve" << std::endl;
////                    for (size_t i = 0; i < count; i++) {
////                        point_t<3> p{ tempData.x(i), tempData.y(i), tempData.z(i) };
////                        points.emplace_back(p);
////                    }
//
//                std::cout << "[" << currentTime() << "] create KDTree start." << std::endl;
//                KDTree<3> tree(tempData);
//                std::cout << "[" << currentTime() << "] create KDTree end." << std::endl;

//                TestData testData;
//                std::cout << "[" << currentTime() << "] create KDTree2 start." << std::endl;
//                KDTree<3> tree2(testData);
//                std::cout << "[" << currentTime() << "] create KDTree2 end." << std::endl;

//                std::cout << "[" << currentTime() << "] test nearest point start." << std::endl;
//                for (int i = 0; i < count; i++) {
//                    auto ret = tree.nearest_index(point_t<3>({ tempData.x(i), tempData.y(i), tempData.z(i) }));
//                }
//                std::cout << "[" << currentTime() << "] test nearest point end." << std::endl;
//
//                std::cout << "[" << currentTime() << "] test neighborhood points start." << std::endl;
//                for (int i = 0; i < count; i++) {
//                    auto ret = tree.neighborhood_indices(point_t<3>({ tempData.x(0), tempData.y(0), tempData.z(0) }), 1);
//                }
//                std::cout << "[" << currentTime() << "] test neighborhood points end." << std::endl;

//                std::cout << "[" << currentTime() << "] test neighborhood points distances start." << std::endl;
//                std::vector<size_t> indexes;
//                std::vector<double> distances;
//                auto ret = tree.neighborhood(point_t<3>({ tempData.x(0), tempData.y(0), tempData.z(0) }), 1, indexes, distances);
//                std::cout << "[" << currentTime() << "] test neighborhood points distances end." << std::to_string(ret) << std::endl;

            } catch (std::exception& exception) {
                std::cout << "err:" << exception.what() << std::endl;
            }
        }
#endif
    }

    static void setLookAtM(std::array<float,16>& rm, float eyeX, float eyeY, float eyeZ,
                           float centerX, float centerY, float centerZ, float upX, float upY,
                           float upZ) {

        // See the OpenGL GLUT documentation for gluLookAt for a description
        // of the algorithm. We implement it in a straightforward way:

        float fx = centerX - eyeX;
        float fy = centerY - eyeY;
        float fz = centerZ - eyeZ;

        // Normalize f
        float rlf = 1.0f / length(fx, fy, fz);
        fx *= rlf;
        fy *= rlf;
        fz *= rlf;

        // compute s = f x up (x means "cross product")
        float sx = fy * upZ - fz * upY;
        float sy = fz * upX - fx * upZ;
        float sz = fx * upY - fy * upX;

        // and normalize s
        float rls = 1.0f / length(sx, sy, sz);
        sx *= rls;
        sy *= rls;
        sz *= rls;

        // compute u = s x f
        float ux = sy * fz - sz * fy;
        float uy = sz * fx - sx * fz;
        float uz = sx * fy - sy * fx;

        rm[0] = sx;
        rm[1] = ux;
        rm[2] = -fx;
        rm[3] = 0.0f;

        rm[4] = sy;
        rm[5] = uy;
        rm[6] = -fy;
        rm[7] = 0.0f;

        rm[8] = sz;
        rm[9] = uz;
        rm[10] = -fz;
        rm[11] = 0.0f;

        rm[12] = 0.0f;
        rm[13] = 0.0f;
        rm[14] = 0.0f;
        rm[15] = 1.0f;

        translateM(rm, -eyeX, -eyeY, -eyeZ);
    }

    inline static float length(float x, float y, float z) {
        return std::sqrt(x * x + y * y + z * z);
    }

    inline static void translateM(
            std::array<float, 16>& m, float x, float y, float z) {
        for (int i=0 ; i<4 ; i++) {
            m[12 + i] += m[i] * x + m[4 + i] * y + m[8 + i] * z;
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
#endif //PLYFASTLOADER_TEST_H
