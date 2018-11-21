#ifndef V8_UTILS_H_
#define V8_UTILS_H_

namespace v8 {
namespace internal {

inline int StrLength(const char* string) {
  size_t length = strlen(string);
  ASSERT(length == static_cast<size_t>(static_cast<int>(length)));
  return static_cast<int>(length);
}

template <typename T>
static T* NewArray(int size) {
  T* result = new T[size];
  assert(result != NULL);
  return result;
}


template <typename T>
static void DeleteArray(T* array) {
  delete[] array;
}

// Returns the maximum of the two parameters.
template <typename T>
static T Max(T a, T b) {
  return a < b ? b : a;
}


// Returns the minimum of the two parameters.
template <typename T>
static T Min(T a, T b) {
  return a < b ? a : b;
}


template <typename T>
static int Compare(const T& a, const T& b) {
  if (a == b)
    return 0;
  else if (a < b)
    return -1;
  else
    return 1;
}


template <typename T>
static int PointerValueCompare(const T* a, const T* b) {
  return Compare<T>(*a, *b);
}


template <class Dest, class Source>
struct BitCastHelper {
  STATIC_ASSERT(sizeof(Dest) == sizeof(Source));

  INLINE(static Dest cast(const Source& source)) {
    Dest dest;
    memcpy(&dest, &source, sizeof(dest));
    return dest;
  }
};

template <class Dest, class Source>
struct BitCastHelper<Dest, Source*> {
  INLINE(static Dest cast(Source* source)) {
    return BitCastHelper<Dest, uintptr_t>::
        cast(reinterpret_cast<uintptr_t>(source));
  }
};

template <class Dest, class Source>
INLINE(Dest BitCast(const Source& source));

template <class Dest, class Source>
inline Dest BitCast(const Source& source) {
  return BitCastHelper<Dest, Source>::cast(source);
}



template <typename T>
class Vector {
 public:
  Vector() : start_(NULL), length_(0) {}
  Vector(T* data, int length) : start_(data), length_(length) {
    ASSERT(length == 0 || (length > 0 && data != NULL));
  }

  static Vector<T> New(int length) {
    return Vector<T>(NewArray<T>(length), length);
  }

  // Returns a vector using the same backing storage as this one,
  // spanning from and including 'from', to but not including 'to'.
  Vector<T> SubVector(int from, int to) {
    ASSERT(to <= length_);
    ASSERT(from < to);
    ASSERT(0 <= from);
    return Vector<T>(start() + from, to - from);
  }

  // Returns the length of the vector.
  int length() const { return length_; }

  // Returns whether or not the vector is empty.
  bool is_empty() const { return length_ == 0; }

  // Returns the pointer to the start of the data in the vector.
  T* start() const { return start_; }

  // Access individual vector elements - checks bounds in debug mode.
  T& operator[](int index) const {
    ASSERT(0 <= index && index < length_);
    return start_[index];
  }

  const T& at(int index) const { return operator[](index); }

  T& first() { return start_[0]; }

  T& last() { return start_[length_ - 1]; }

  // Returns a clone of this vector with a new backing store.
  Vector<T> Clone() const {
    T* result = NewArray<T>(length_);
    for (int i = 0; i < length_; i++) result[i] = start_[i];
    return Vector<T>(result, length_);
  }

  void Sort(int (*cmp)(const T*, const T*)) {
    typedef int (*RawComparer)(const void*, const void*);
    qsort(start(),
          length(),
          sizeof(T),
          reinterpret_cast<RawComparer>(cmp));
  }

  void Sort() {
    Sort(PointerValueCompare<T>);
  }

  void Truncate(int length) {
    ASSERT(length <= length_);
    length_ = length;
  }

  // Releases the array underlying this vector. Once disposed the
  // vector is empty.
  void Dispose() {
    DeleteArray(start_);
    start_ = NULL;
    length_ = 0;
  }

  inline Vector<T> operator+(int offset) {
    ASSERT(offset < length_);
    return Vector<T>(start_ + offset, length_ - offset);
  }

  // Factory method for creating empty vectors.
  static Vector<T> empty() { return Vector<T>(NULL, 0); }

  template<typename S>
  static Vector<T> cast(Vector<S> input) {
    return Vector<T>(reinterpret_cast<T*>(input.start()),
                     input.length() * sizeof(S) / sizeof(T));
  }

 protected:
  void set_start(T* start) { start_ = start; }

 private:
  T* start_;
  int length_;
};

static inline char* StrDup(const char* str) {
  int length = StrLength(str);
  char* result = NewArray<char>(length + 1);
  memcpy(result, str, length);
  result[length] = '\0';
  return result;
}

} }  // namespace v8::internal

#endif  // V8_UTILS_H_
