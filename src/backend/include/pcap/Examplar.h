#include <string>

using namespace std;

class Examplar {
private:
	int _size;
	float* _map;
public:
	Examplar() {
		_map = NULL;
	}
		
	~Examplar() {
		printf("deleting\n");
		if (_map) {
			printf("real deleting\n");
			delete[] _map;
		}
	}

	void SetBuffer(int size) {
		_size = size;
		_map = new float[size*size];
	}	
	
	inline float GetI(int x, int y) const {
		return _map[x*_size+y];
	}
	inline void SetI(int x, int y, float imgI) {
		_map[x*_size+y] = imgI;
	}
	inline int size() {
		return _size;
	}
};
