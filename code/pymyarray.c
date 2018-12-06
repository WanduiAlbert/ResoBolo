#include <Python.h>
#include "mylib.c"

/* Define the PyMyarray object structure */
typedef struct{
    PyObject_HEAD;
    /* Type-specific fields go below */
    MyArr arr;
} PyMyArray;

/* __init__ function implemented in C */
static int PyMyArray_init(PyMyArray *self, PyObject *args, PyObject *kwds)
{
    // check if init has already been called
    if (self -> arr.arr != NULL;
        deallocate_MyArray(&self -> arr);

    int length = 0;
    static char *kwlist[] = {"length", NULL};
    if (! PyArgs_ParseTupleAndKeywords(args, kwds, "|i", kwlist, &length))
        return -1;

    if (length < 0)
        length = 0;

    initialize_MyArray(&self -> arr, length);

    return 0;

}

/* function called when the array is deallocated */
static void PyMyArray_dealloc(PyMyArray *self)
{
    deallocate_MyArray(&self -> arr);
    Py_TYPE(self) -> tp_free((PyObject *)self);
}


