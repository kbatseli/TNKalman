A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification (Matlab&copy;/Octave&copy;)
-------------------------------------------------------------------------------------------------------------------------------

The Tensor Network Kalman filter can estimate state vectors that are exponentially large without ever having to explicitly construct them. It also easily accommodates the case where several different state vectors need to be estimated simultaneously. The standard Kalman equations are rewritten as tensor equations and then implemented using Tensor Networks, which effectively transforms the exponential storage cost and computational complexity into a linear one. Once can apply this TN Kalmam filter for recursive nonlinear system identification of high-order discrete-time multiple-input multiple-output (MIMO) Volterra systems.

Tensor Networks are implemented as structures with two fields: *core* and *n*. The *core* field contains a cell of TN-cores. The *n* field contains a matrix of TN-ranks, where size(n,1) is the number of cores and size(n,2) indicates the order of the TN-core. In order to maintain consistency the values of n(1,1) and n(end,end) always need to be 1.

1. Functions
------------

* [m,P]=TNKalman(m,P,A,Q,r,y,u,tol)

Computes the predict and update step for both the mean vectors m and covariance matrices P in the Tensor Network format specifically for the recursive identification of Volterra systems.

* c=addTN(a,b)

Adds two Tensor Networks *a* and *b* together.

* a=cmodeprod(a,u,k,l)

Computes the k-mode product of the lth core of a Tensor Network *a* with the matrix *u*.

* b=contract(a)

Sums the Tensor Network *a* over all its auxiliary indices to obtain the underlying tensor.

* s=frobnorm(a)

Computes the Frobenius norm of a Tensor Network *a*.

* m=initm(l,n,d)

Constructs a Tensor Network *m* corresponding with a n^d x l zero matrix. Use this function to initialize the matrix of means M(0) of the linear state model.

* P=initP(sigmas,n,d)

Constructs a Tensor Network *P* corresponding with a n^d x n^d x l tensor of diagonal matrix slices. Use this function to initialize the tensor of diagonal covariance matrices P(0) of the linear state model.

* a=modeprod(a,u,k)

Computes the k-mode product of all cores of the Tensor Network *a* with the matrix *u*.

* C=khatri(A,B)

Naive implementation of the column-wise Kronecker (Khatri-Rao) product of two matrices *A* and *B*.

* C=mkhatri(A,d)

Calls the khatri.m function d-times on the matrix A.

* a=roundTN(a,tol)

Returns an approximation of the Tensor Network *a* such that the approximation has a relative error *tol*.

* yhat=sim_volterraTN(u,TN)

Simulates a MIMO Volterra system *TN*  in the Tensor Network format with a given input *u*.

* c=squareTN(a,b)

Computes the Tensor Network $c$ of the column-wise outer product of the two matrices $a$ and $b$ in the Tensor Network format.

* a=squeeze(a)

Removes redundant singleton dimensions from a Tensor Network *a*.

* demo.m

Small demo that illustrates how to use the TN Kalman filter for recursive system identifcation of Volterra systems.


2. Reference
------------

"A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification"

Authors: Kim Batselier, Ngai Wong
