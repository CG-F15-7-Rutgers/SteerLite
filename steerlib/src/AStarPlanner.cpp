//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//


#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm> 
#include <functional>
#include <queue>
#include <math.h>
#include "planning/AStarPlanner.h"


#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define EUCLEDIAN 1 // 0 for Manhattan distance or 1 for Eucledian
#define WEIGHT 8 // weight used in heuristic fuction

namespace SteerLib
{
	AStarPlanner::AStarPlanner(){}

	AStarPlanner::~AStarPlanner(){}

	bool AStarPlanner::canBeTraversed ( int id ) 
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x,z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		x_range_min = MAX(x-OBSTACLE_CLEARANCE, 0);
		x_range_max = MIN(x+OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsX());

		z_range_min = MAX(z-OBSTACLE_CLEARANCE, 0);
		z_range_max = MIN(z+OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsZ());


		for (int i = x_range_min; i<=x_range_max; i+=GRID_STEP)
		{
			for (int j = z_range_min; j<=z_range_max; j+=GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords( i, j );
				traversal_cost += gSpatialDatabase->getTraversalCost ( index );
				
			}
		}

		if ( traversal_cost > COLLISION_COST ) 
			return false;
		return true;
	}



	Util::Point AStarPlanner::getPointFromGridIndex(int id)
	{
		Util::Point p;
		gSpatialDatabase->getLocationFromIndex(id, p);
		return p;
	}

	double estimateCost(AStarPlannerNode agentStart, AStarPlannerNode agentGoal)
	{
		double cost;
		if (EUCLEDIAN)
		{
			cost = sqrt(pow((agentStart.point.x-agentGoal.point.x),2)+pow((agentStart.point.z-agentGoal.point.z),2));
		}
		else
		{
			cost = abs(agentStart.point.x - agentGoal.point.x) + abs(agentStart.point.z - agentGoal.point.z);
		}
		return cost;
	}
	
	AStarPlannerNode findLowFScore(std::vector<AStarPlannerNode> &open, int &index)
	{
		index = INT_MAX;
		AStarPlannerNode node;
		std::vector<AStarPlannerNode>::iterator itr = open.begin();
		for (int i=0; itr != open.end(); itr++,i++)
		{
			if (node > *itr) {
				index = i;
				node = *itr;
			}
			if (node.f == (*itr).f)
			{
				if (node.g < (*itr).g)
				{
					node = *itr;
					index = i;
				}
			}

		}
		return node;
	}

	std::vector<AStarPlannerNode> AStarPlanner::getNeighbours(AStarPlannerNode currentNode, AStarPlannerNode goalNode)
	{
		int currentIndex = gSpatialDatabase->getCellIndexFromLocation(currentNode.point);
		std::vector<AStarPlannerNode> neighbours;
		unsigned int currentPointx,currentPointz; 
		gSpatialDatabase->getGridCoordinatesFromIndex(currentIndex, currentPointx, currentPointz);
		Util::Point currentPoint(currentPointx, 0, currentPointz);
		Util::Point nextPoint(0,0,0);
		for (int i = 0; i < 8; i++)
		{
			nextPoint = currentPoint;

			switch (i)
			{
				case 0:
					nextPoint.x++;
					break;
				case 1:
					nextPoint.x--;
					break;
				case 2:
					nextPoint.z++;
					break;
				case 3:
					nextPoint.z--;
					break;
				case 4:
					nextPoint.x++;
					nextPoint.z++;
					break;
				case 5:
					nextPoint.x++;
					nextPoint.z--;
					break;
				case 6:
					nextPoint.x--;
					nextPoint.z--;
					break;
				case 7:
					nextPoint.x--;
					nextPoint.z++;
					break;
			}

			int index = gSpatialDatabase->getCellIndexFromGridCoords(nextPoint.x, nextPoint.z);
			if (canBeTraversed(index))
			{
				AStarPlannerNode node;
				nextPoint = getPointFromGridIndex(index);
				node.point = nextPoint;
				node.g = currentNode.g + 1;
				node.f = node.g + (WEIGHT*estimateCost(node, goalNode));
				node.parentx = currentNode.point.x;
				node.parentz = currentNode.point.z;
				neighbours.push_back(node);
			}
		
		}
	
		return neighbours;
	}

	bool closedList(AStarPlannerNode nodeToCheck, std::vector<AStarPlannerNode> closedNodes)
	{
		for (std::vector<AStarPlannerNode>::iterator itr = closedNodes.begin(); itr != closedNodes.end(); itr++)
		{
			if (nodeToCheck == (*itr))
				return true;
		}
		return false;
	}

	int openList(AStarPlannerNode nodeToCheck, std::vector<AStarPlannerNode> closedNodes)
	{
		int i = 0;
		for (std::vector<AStarPlannerNode>::iterator itr = closedNodes.begin(); itr != closedNodes.end(); itr++,i++)
		{
			if (nodeToCheck == (*itr))
				return i;
		}
		i = -1;
		return i;
	}

	AStarPlannerNode findClosedNodes(std::vector<AStarPlannerNode> closedNodes, AStarPlannerNode nodeToFind)
	{
		for (std::vector<AStarPlannerNode>::iterator itr = closedNodes.begin(); itr != closedNodes.end(); itr++)
		{
			if ((*itr).point.x == nodeToFind.parentx && (*itr).point.z == nodeToFind.parentz)
				return (*itr);
		}
	}

	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path,  Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path)
	{
		gSpatialDatabase = _gSpatialDatabase;

		int startIndex =_gSpatialDatabase->getCellIndexFromLocation(start);
		start = getPointFromGridIndex(startIndex);
		AStarPlannerNode Astart(start, 0, 0, NULL);
		int goalIndex = _gSpatialDatabase->getCellIndexFromLocation(goal);
		goal = getPointFromGridIndex(goalIndex);
		AStarPlannerNode Agoal(goal,DBL_MAX,DBL_MAX,NULL);
		std::vector<AStarPlannerNode> open;
		std::vector<AStarPlannerNode> closed;
		Astart.g = 0;
		Astart.f = Astart.g + estimateCost(Astart,Agoal);
		open.emplace_back(Astart);
		int index;
		std::vector<AStarPlannerNode> neighbours;
		AStarPlannerNode current;
		while (!open.empty())
		{
			neighbours.clear();
		current = findLowFScore(open,index);
		if (current.point == Agoal.point)
		{   			
			std::cout << "\nPath length:" << current.g;

			while (!(current.parentx == 0 && current.parentz==0))
			{
				agent_path.push_back(current.point);
				current = findClosedNodes(closed, current);
			}
			std::reverse(agent_path.begin(),agent_path.end());
			std::cout << "\n\tNumber of nodes in Closed List: " << closed.size();
			return true;
		}
		open.erase(open.begin()+index);
		
		closed.push_back(current);

		neighbours = getNeighbours(current,Agoal);
		double tentative_gscore=0;
		int i; int z = 0;
		for (std::vector<AStarPlannerNode>::iterator itr = neighbours.begin(); itr != neighbours.end(); itr++,z++)
		{
			if (closedList(*itr, closed))
				continue;
			tentative_gscore = current.g + 1;
			i = openList(*itr, open);
			if (i != -1)
			{
				if (tentative_gscore < open.at(i).g)
				{   
					open.at(i).g = tentative_gscore;
					open.at(i).parentx = current.point.x;
					open.at(i).parentz = current.point.z;
					open.at(i).f = open.at(i).g + estimateCost(open.at(i),Agoal);
				}
			}
			else
			{
				open.push_back(neighbours.at(z));
			}
		}
		
		}
		std::cout << "\n No path found";

		return false;
	}
}
