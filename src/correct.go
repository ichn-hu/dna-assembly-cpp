package main

func readFasta(path string) -> string[] {
	
}

func main() {
    s := "Hello Cgo"
    cs := C.CString(s)
    C.print(cs)
    C.free(unsafe.Pointer(cs))
}