/*!
 *
 * \author VaHiD AzIzI
 *
 */


#include "obstacles/GJK_EPA.h"

// not needed because we just call gjk
SteerLib::GJK_EPA::GJK_EPA()
{
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
    std::vector<Util::Vector> simplex;
    
    bool isColliding = GJK(_shapeA, _shapeB, simplex);
    
    if (isColliding) {
        epa(_shapeA, _shapeB, simplex, return_penetration_depth, return_penetration_vector);
    }
    
    return isColliding;
}

bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB, std::vector<Util::Vector>& simplex) {

    Util::Vector centerShapeA(0,0,0);
    Util::Vector centerShapeB(0,0,0);
    Util::Vector d(0,0,0);

    centerShapeA = central_origin_of_shape(shapeA);
    centerShapeB = central_origin_of_shape(shapeB);
   
    // direction
    d = centerShapeB - centerShapeA;
    
    simplex.push_back(support_function(shapeA, shapeB, d));
    
    // negate the direction
    d = -d;
    
    while (true) {
        
        simplex.push_back(support_function(shapeA, shapeB, d));
        
        if (dot(simplex.back(), d) <= 0){
            return false;
        }
        else if (containsOrigin(simplex, d)) {
            
            if (simplex.size() < 3) {
                simplex.push_back(support_function(shapeA, shapeB, d));
            }
            
            return true;
        }
    }
}

bool SteerLib::GJK_EPA::containsOrigin(std::vector<Util::Vector>& simplex, Util::Vector& direction) {
 // have to declare loads of variables for this function!!
    Util::Vector first_point_in_simplex(0,0,0);
    Util::Vector negation_of_1st_point (0,0,0);
    Util::Vector pointB (0,0,0);
    Util::Vector pointC (0,0,0);
    Util::Vector abEdge (0,0,0);
    Util::Vector acEdge (0,0,0);
    Util::Vector abPerp (0,0,0);
    Util::Vector acPerp (0,0,0);
    
    first_point_in_simplex = simplex.back();
    negation_of_1st_point = -first_point_in_simplex;
    
    // this means it is a triangle
    if (simplex.size() == 3) {
        
         pointB = simplex[1];
         pointC = simplex[0];
        
        abEdge = pointB - first_point_in_simplex;
        acEdge = pointC - first_point_in_simplex;
        
        abPerp = abEdge * dot(abEdge, acEdge) - acEdge * (dot(abEdge, abEdge));
        acPerp = acEdge * dot(acEdge, abEdge) - abEdge * (dot(acEdge, acEdge));
     
        if (dot(abPerp, negation_of_1st_point) > 0) {
            simplex.erase(simplex.begin());
            direction = abPerp;
        }
        else if (dot(acPerp, negation_of_1st_point) > 0) {
            
            simplex.erase(simplex.begin() + 1);
            direction = acPerp;
        }
        else return true;
    }
    else {
        
        Util::Vector pointB = simplex.at(0);
        Util::Vector abEdge = pointB - first_point_in_simplex;
        Util::Vector abPerp = negation_of_1st_point * dot(abEdge, abEdge) - abEdge * dot(abEdge, negation_of_1st_point);
        direction = abPerp;
        if (dot(abPerp, negation_of_1st_point) == 0) {
            float aToBDotnegation_of_1st_point = dot(abEdge, negation_of_1st_point);
            if (aToBDotnegation_of_1st_point >= 0 && aToBDotnegation_of_1st_point < dot(abEdge, abEdge)) return true;
        }
    }
    return false;
}

Util::Vector SteerLib::GJK_EPA::Minkowski_Helper(const std::vector<Util::Vector>& shape, const Util::Vector& direction) {
    
    Util::Vector furthestPoint(0, 0, 0);
    float farthestDistance = dot(shape[0], direction);
    int farthestIndex = 0;

    int i = 1;
    while (i < shape.size()) {
        float dotProduct = dot(shape[i], direction);
        if (dotProduct > farthestDistance) {
            farthestDistance = dotProduct;
            farthestIndex = i;
        }
        i++;
    }
    
    furthestPoint[0] = shape[farthestIndex][0];
    furthestPoint[1] = shape[farthestIndex][1];
    furthestPoint[2] = shape[farthestIndex][2];
    
    return furthestPoint;
}

// support function, returns a vector of the minkowski difference
Util::Vector SteerLib::GJK_EPA::support_function(const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB, const Util::Vector& direction) {
    Util::Vector Minkowski_Difference = Minkowski_Helper(shapeA, direction) - Minkowski_Helper(shapeB, -direction);
    return Minkowski_Difference;
}

Util::Vector SteerLib::GJK_EPA::central_origin_of_shape(const std::vector<Util::Vector>& shape) {
    
    Util::Vector the_shapes_center (0, 0, 0);
    Util::Vector central_location  (0, 0, 0);
    Util::Vector temporary_vector  (0, 0, 0);
    
    float shapeSize = shape.size();
    int count = 0;
    
    float xvalue = 0;
    float zvalue = 0;
    
    while(count < shapeSize)
    {
        central_location = shape[count];
        
        xvalue = central_location[0];
        zvalue = central_location[2];
        
        the_shapes_center [0] = the_shapes_center [0] + xvalue;
        the_shapes_center [2] = the_shapes_center [2] + xvalue;
        
       temporary_vector [0] = the_shapes_center [0] / shapeSize;
       temporary_vector [2] = the_shapes_center [2] / shapeSize;
        
     
        count++;
    }
    
    the_shapes_center = temporary_vector;
    
    return the_shapes_center;
}

float SteerLib::GJK_EPA::dot(const Util::Vector& vectorA, const Util::Vector& vectorB) {
    float dot_product = (vectorA[0]* vectorB[0])+(vectorA[1]* vectorB[1])+(vectorA[2]* vectorB[2]);
    return dot_product;
}

void SteerLib::GJK_EPA::epa(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex, float& penetration_depth, Util::Vector& penetration_vector)
{
    Util::Vector normal (0,0,0);
    Util::Vector supportPoint (0,0,0);
    // has to be a point close to 0, but not actually zero.
    float TOLERANCE = 0.000001;
    float distance = 0;
    int i  = 0;
    
    while (true)
    {
        findNearestEdge(simplex, distance, normal, i);
        
        supportPoint = support_function(_shapeA, _shapeB, normal);
        
        double d = supportPoint * normal;
        
        if (d - distance < TOLERANCE) {
            penetration_vector = normal;
            penetration_depth = d;
            break;
            
        } else {
            simplex.insert(simplex.begin()+i, supportPoint);
        }
    }
    
    return;
}

void SteerLib::GJK_EPA::findNearestEdge(std::vector<Util::Vector> simplex, float& distance, Util::Vector& normal, int& index)
{
    Util::Vector current_point;
    Util::Vector point_plus_one;
    Util::Vector e;
    Util::Vector n;
    Util::Vector n_norm;

    distance = FLT_MAX;

    int count;
    for (int count = 0; count < simplex.size(); count++) {
        int j = count + 1 == simplex.size() ? 0 : count + 1;
        current_point = simplex[count];
        point_plus_one = simplex[j];

        e = point_plus_one - current_point;
        n = tripleProduct(e, current_point, e);
        n_norm = n / n.norm();

        double d = n_norm * current_point;

        if (d < distance) {
            distance = d;
            normal = n_norm;
            index = j;
        }
    }
}

Util::Vector SteerLib::GJK_EPA::tripleProduct(Util::Vector A, Util::Vector B, Util::Vector C)
{
    return (B * (C * A)) - (A * (C * B));
}
