#include <string>

using namespace std;

class Center {
private:
	int x, y, radius, im_id;
public:
	Center();
	Center(int x, int y) {
		this->x = x;
		this->y = y;
	}
	
	inline int GetX() const {
		return x;
	}
	inline int GetY() const {
		return y;
	}
	inline int GetR() const {
		return radius;
	}
	inline int GetID() const {
		return im_id;
	}
	inline void SetX(int x) {
		this->x = x;
	}
	inline void SetY(int y) {
		this->y = y;
	}
	inline void SetR(int r) {
		this->radius = r;
	}
	inline void SetID(int id) {
		this->im_id = id;
	}
};
