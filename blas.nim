when defined(windows):
  const blasSuffix = ".dll"
elif defined(macosx):
  const blasSuffix = ".dylib"
else:
  const blasSuffix = ".so"

const blasPrefix = "libblas"
const blasName = blasPrefix & blasSuffix

import acomplex

proc gemm*(
  transa, transb: ptr char,
  m,n,k: ptr int,
  alpha,a: ptr float32,
  lda: ptr int,
  b: ptr float64,
  ldb: ptr int,
  beta,c: ptr float32,
  ldc: ptr int,
) {.cdecl, importc: "sgemm_", dynlib: blasName.}

proc gemm*(
  transa, transb: ptr char,
  m,n,k: ptr int,
  alpha,a: ptr float64,
  lda: ptr int,
  b: ptr float64,
  ldb: ptr int,
  beta,c: ptr float64,
  ldc: ptr int,
) {.cdecl, importc: "dgemm_", dynlib: blasName.}

proc gemm*(
  transa, transb: ptr char,
  m,n,k: ptr int,
  alpha,a: ptr complex32,
  lda: ptr int,
  b: ptr complex32,
  ldb: ptr int,
  beta,c: ptr complex32,
  ldc: ptr int,
) {.cdecl, importc: "cgemm_", dynlib: blasName.}

proc gemm*(
  transa, transb: ptr char,
  m,n,k: ptr int,
  alpha,a: ptr complex64,
  lda: ptr int,
  b: ptr complex64,
  ldb: ptr int,
  beta,c: ptr complex64,
  ldc: ptr int,
) {.cdecl, importc: "zgemm_", dynlib: blasName.}
