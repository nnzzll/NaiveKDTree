#pragma once
#include <algorithm>
#include <iostream>
#include <set>
#include <stack>
#include <vector>

struct KDNode
{
	int index;
	int axis;
	KDNode *left;
	KDNode *right;
	KDNode(int index, int axis, KDNode *left = nullptr, KDNode *right = nullptr)
	{
		this->index = index;
		this->axis = axis;
		this->left = left;
		this->right = right;
	}
};

template <class T>
class KDTree
{
private:
	int ndim;
	KDNode *root;
	KDNode *build(std::vector<T> &);
	std::set<int> visited;
	std::stack<KDNode *> queueNode;
	std::vector<T> m_data;

	void release(KDNode *);
	void printNode(KDNode *);
	int chooseAxis(std::vector<T> &);
	void dfs(KDNode *, T);
	// 点与点之间的距离
	inline double distanceT(KDNode *, T);
	inline double distanceT(int, T);
	// 点与超平面的距离
	inline double distanceP(KDNode *, T);
	// 检查父节点超平面是否在超球体中
	inline bool checkParent(KDNode *, T, double);

public:
	KDTree(std::vector<T> &, int);
	~KDTree();
	void Print();
	int findNearestPoint(T);
};

template <class T>
KDTree<T>::KDTree(std::vector<T> &data, int dim)
{
	ndim = dim;
	m_data = data;
	root = build(data);
}

template <class T>
KDTree<T>::~KDTree()
{
	release(root);
}

template <class T>
void KDTree<T>::Print()
{
	printNode(root);
}

template <class T>
void KDTree<T>::release(KDNode *node)
{
	if (node)
	{
		if (node->left)
			release(node->left);
		if (node->right)
			release(node->right);
		delete node;
		node = nullptr;
	}
}

template <class T>
KDNode *KDTree<T>::build(std::vector<T> &data)
{
	if (data.empty())
		return nullptr;
	std::vector<T> temp = data;
	int mid_index = static_cast<int>(data.size() / 2);
	int axis = data.size() > 1 ? chooseAxis(temp) : -1;
	std::sort(temp.begin(), temp.end(), [axis](T a, T b)
			  { return a[axis] < b[axis]; });
	std::vector<T> leftData, rightData;
	leftData.assign(temp.begin(), temp.begin() + mid_index);
	rightData.assign(temp.begin() + mid_index + 1, temp.end());
	KDNode *leftNode = build(leftData);
	KDNode *rightNode = build(rightData);
	KDNode *rootNode;
	rootNode = new KDNode(temp[mid_index].index, axis, leftNode, rightNode);
	return rootNode;
}

template <class T>
void KDTree<T>::printNode(KDNode *node)
{
	if (node)
	{
		std::cout << "Index: " << node->index << "\tAxis: " << node->axis << std::endl;
		printNode(node->left);
		printNode(node->right);
	}
}

template <class T>
int KDTree<T>::chooseAxis(std::vector<T> &data)
{
	int axis = -1;
	double maxVar = -1.0;
	size_t size = data.size();
	for (int i = 0; i < ndim; i++)
	{
		double mean = 0;
		double Var = 0;
		for (int j = 0; j < size; j++)
		{
			mean += static_cast<double>(data[j][i]);
		}
		mean = mean / static_cast<double>(size);
		for (int j = 0; j < size; j++)
		{
			Var += (static_cast<double>(data[j][i]) - mean) * (static_cast<double>(data[j][i]) - mean);
		}
		Var = Var / static_cast<double>(size);
		if (Var > maxVar)
		{
			axis = i;
			maxVar = Var;
		}
	}
	return axis;
}

template <class T>
int KDTree<T>::findNearestPoint(T pt)
{
	while (!queueNode.empty())
		queueNode.pop();
	double min_dist = DBL_MAX;
	int resNodeIdx = -1;
	dfs(root, pt);
	while (!queueNode.empty())
	{
		KDNode *curNode = queueNode.top();
		queueNode.pop();
		double dist = distanceT(curNode, pt);
		if (dist < min_dist)
		{
			min_dist = dist;
			resNodeIdx = curNode->index;
		}

		if (!queueNode.empty())
		{
			KDNode *parentNode = queueNode.top();
			int parentAxis = parentNode->axis;
			int parentIndex = parentNode->index;
			if (checkParent(parentNode, pt, min_dist))
			{
				if (m_data[curNode->index][parentNode->axis] < m_data[parentNode->index][parentNode->axis])
					dfs(parentNode->right, pt);
				else
					dfs(parentNode->left, pt);
			}
		}
	}
	return resNodeIdx;
}

template <class T>
void KDTree<T>::dfs(KDNode *node, T pt)
{
	if (node)
	{
		if (visited.find(node->index) != visited.end())
			return;
		queueNode.push(node);
		visited.insert(node->index);
		if (pt[node->axis] < m_data[node->index][node->axis] && node->left)
			return dfs(node->left, pt);
		else if (pt[node->axis] >= m_data[node->index][node->axis] && node->right)
			return dfs(node->right, pt);
	}
}

template <class T>
double KDTree<T>::distanceT(KDNode *node, T pt)
{
	double dist = 0;
	for (int i = 0; i < ndim; i++)
	{
		dist += (m_data[node->index][i] - pt[i]) * (m_data[node->index][i] - pt[i]);
	}
	return sqrt(dist);
}

template <class T>
double KDTree<T>::distanceT(int index, T pt)
{
	double dist = 0;
	for (int i = 0; i < ndim; i++)
	{
		dist += (m_data[index][i] - pt[i]) * (m_data[index][i] - pt[i]);
	}
	return sqrt(dist);
}

template <class T>
double KDTree<T>::distanceP(KDNode *node, T pt)
{
	int axis = node->axis;
	double dist = static_cast<double>(pt[axis] - m_data[node->index][axis]);
	return dist >= 0 ? dist : -dist;
}

template <class T>
bool KDTree<T>::checkParent(KDNode *node, T pt, double distT)
{
	double distP = distanceP(node, pt);
	return distP <= distT;
}