#ifndef FRAMEWORK_COMMON_MACRO_H__
#define FRAMEWORK_COMMON_MACRO_H__

// disable unuse warning
#define DD_UNUSED(x) (void)x

/*! concatenating multiple args into one*/
#define CONCATE(...) __VA_ARGS__

/*! getter and setter generator for class memeber */
#define ADD_CLASS_FIELD(type, name, getter, setter) \
	public: \
	type& getter() { return _##name; } \
	type const & getter() const{ return _##name; } \
	void setter(const type& name) { _##name = name; } \
	private: \
	type _##name;

#define ADD_CLASS_FIELD_PRIVATE_SETTER(type, name, getter, setter) \
	public: \
	type const & getter() const{ return _##name; } \
	private: \
	void setter(const type& name) { _##name = name; } \
	private: \
	type _##name;

#define ADD_CLASS_FIELD_NOSETTER(type, name, getter) \
	public: \
	type& getter() { return _##name; } \
	type const & getter() const{ return _##name; } \
	private: \
	type _##name;

#define ADD_CLASS_FIELD_NOGETTER(type, name, setter) \
	public: \
	void setter(const type& name) { _##name = name; } \
	private: \
	type _##name;

#define ADD_CLASS_FIELD_PRIVATE(type, name ) \
	private: \
	type _##name;

#define ADD_POINTER_FIELD_NOSET(type, name, getter) \
	public: \
	type& getter() { return *_##name; } \
	const type & getter() const{ return *_##name; } \
	private: \
	type* _##name;


#define SAFE_NEW_POINTER(p, T, ...)         \
if (!p)                                 \
	{                                       \
	p = new T(__VA_ARGS__);             \
	}                                       \

#define SAFE_DELETE_POINTER(p)              \
if (p)                                  \
	{                                       \
	delete p;                           \
	p = NULL;                           \
	}


// two uint8_t transfer with unsignd short
#define FW_USHORT_MAKE(high, low)	((uint16_t)(((uint8_t)((low) & 0xFF)) | ((uint16_t)((uint8_t)((high) & 0xFF))) << 8))
#define FW_USHORT_LOW(i)            ((uint8_t)((uint16_t)(i) & 0xFF))
#define FW_USHORT_HIGH(i)           ((uint8_t)((uint16_t)(i) >> 8))

// two uint16_t transfer with uint32_t
#define FW_UINT_MAKE(high, low)		((uint32_t)(((uint16_t)((low) & 0xFFFF)) | ((uint32_t)((uint16_t)((high) & 0xFFFF))) << 16))
#define FW_UINT_LOW(i)				((uint16_t)((uint32_t)(i) & 0xFFFF))
#define FW_UINT_HIGH(i)				((uint16_t)((uint32_t)(i) >> 16))

// two uint32_t transfer with uint64
#define FW_UINT64_MAKE(high, low)   ((uint64_t)(((uint32_t)((low) & 0xFFFFFFFF)) | ((uint64_t)((uint32_t)((high) & 0xFFFFFFFF))) << 32))
#define FW_UINT64_LOW(i)            ((uint32_t)((uint64_t)(i) & 0xFFFFFFFF))
#define FW_UINT64_HIGH(i)           ((uint32_t)((uint64_t)(i) >> 32))

// 32bit  bit operation
// pos' value must be 0-31
#define FW_BIT_SET(i, pos)			(i |= ((uint32_t)1 << pos))        // set the pos bit by 1
#define FW_BIT_CLR(i, pos)			(i &= ~((uint32_t)1 << pos))       // set the pos bit by 0
#define FW_BIT_TEST(i, pos)			((i & ((uint32_t)1 << pos)) != 0)  // check the pos bit is 1

#endif  // FRAMEWORK_COMMON_MACRO_H__

