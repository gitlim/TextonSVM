#include <vector>

using namespace std;

template <typename T> class safe_vector {
private:
	vector<T> vec;
public:
	inline int size() const {
		return vec.size();
	}	
	inline void push_back(T e) {
		vec.push_back(e);
	}
	inline void pop_back() {
		vec.pop_back();
	}
	inline T back() const {
		return vec.back();
	}
	inline T at(int i) const {
		return vec[i];
	}
	inline void reset() {
		while (!vec.empty()) {
			delete vec.back();
			vec.pop_back();
		}
	}
	~safe_vector() {
		reset();
	}
};
