import acomplex
import lapack
import blas

type
  Matrix*[M, N: static[int], T] = array[M, array[N, T]]

proc `$`*[M,N:static[int],T](a: Matrix[M,N,T]): string =
  result = ""
  for m in 0..<M:
    if m>0:
      result &= "\n"
    for n in 0..<N:
      result &= $a[m][n]&' '

proc transpose*[M,N:static[int],T](a: Matrix[M,N,T]): Matrix[N,M,T] =
  for m in 0..<M:
    for n in 0..<N:
      result[n][m] = a[m][n]

proc eye*[T](N: static[int], d:T): auto =
  var res: array[N, array[N, T]]
  # var res: Matrix[N,N,T]
  for i in 0..<N:
    for j in 0..<N:
      if i == j:
        res[i][j] = T(d)
      else:
        res[i][j] = T(0)
  return res

proc fill*[T](M: static[int], N: static[int], d:T): auto =
  var res: array[M, array[N, T]]
  # var res: Matrix[M,N,T]
  for i in 0..<M:
    for j in 0..<N:
      res[i][j] = T(d)
  return res

proc `+`*[M,N:static[int],T](a,b: Matrix[M,N,T]): Matrix[M,N,T] =
  for m in 0..<M:
    for n in 0..<N:
      result[m][n] = a[m][n]+b[m][n]

proc `*`*[M,N:static[int],T](a: Matrix[M,N,T], b: SomeNumber): Matrix[M,N,T] =
  for m in 0..<M:
    for n in 0..<N:
      result[m][n] = a[m][n]*T(b)

proc `*`*[M,N:static[int],T](a: SomeNumber, b: Matrix[M,N,T]): Matrix[M,N,T] =
  for m in 0..<M:
    for n in 0..<N:
      result[m][n] = T(a)*b[m][n]

proc `-`*[M,N:static[int],T](a: Matrix[M,N,T]): Matrix[M,N,T] =
  -1*a

proc `-`*[M,N:static[int],T](a,b: Matrix[M,N,T]): Matrix[M,N,T] =
  a+(-b)

proc `*`*[M,N,K:static[int],T](a: var Matrix[M,K,T], b: var Matrix[K,N,T]): Matrix[M,N,T] =
  var
    transa = 'N'
    transb = 'N'
    m = M.int
    n = N.int
    k = K.int
  when T is complex:
    var alpha = [1.0, 0.0]
    var beta = [0.0, 0.0]
  else:
    var alpha = 1.0
    var beta = 0.0
  gemm(transb.addr, transa.addr, n.addr, m.addr, k.addr, alpha.addr, b[0][0].addr, n.addr, a[0][0].addr, k.addr, beta.addr, result[0][0].addr, n.addr)

proc `\`*[M,N:static[int],T](a: var Matrix[M,M,T], b: var Matrix[M,N,T]): Matrix[M,N,T] =
  var
    m = M.int
    n = N.int
    ta = a.transpose
    tb = b.transpose
    ipiv: array[M, int]
    nfo: int
  gesv(m.addr, n.addr, ta[0][0].addr, m.addr, ipiv[0].addr, tb[0][0].addr, m.addr, nfo.addr)
  if nfo != 0:
    info(nfo, "gesv")
  else:
    result = tb.transpose

proc `inv`*[N:static[int],T](a: var Matrix[N,N,T]): Matrix[N,N,T] =
  var
    tb: Matrix[N,N,T] = eye(N,T(1))
    m = N.int
    n = N.int
    ta = a.transpose
    ipiv: array[N, int]
    nfo: int
  gesv(m.addr, n.addr, ta[0][0].addr, m.addr, ipiv[0].addr, tb[0][0].addr, m.addr, nfo.addr)
  if nfo != 0:
    info(nfo, "gesv")
  else:
    result = tb.transpose
