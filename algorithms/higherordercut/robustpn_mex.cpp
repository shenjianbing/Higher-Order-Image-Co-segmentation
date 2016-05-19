#include <mex.h>
#include <stdio.h>
#include <string.h>
#include "jm.h"
#include "weight.h"
#include "expand.h" // minimize energy using alpha-expantion steps

// declarations
template<typename termType>
uint *graphchop(uint dims[3], uint N, uint label_cnt,
                t_weight_fn R_fn, n_weight_fn B_fn, double *V, Energy<termType> *e);
template<class T>
void GetArr(const mxArray* x, T* arr, T bias = 0);

// consts
const int HOP_N_OF_FIELDS(4); // expecting 4 fields for the HOpotentials struct
const char* HOP_FIELDS[HOP_N_OF_FIELDS] = {"ind", "w", "gamma", "Q"};
const int MAX_ITER(50); // maximum 50 iterations
#define sub2ind(r,c,d) (dims[0]*dims[1]*(d) + dims[0]*(c) + (r))

/* pre-calculated edge families */

/*-- 2D --*/
int N4[13][3] = {{ 0, 1, 0 },
              { -1,0, 0 }};
int N8[13][3] = {{ 0, 1, 0 },
              { -1, 1, 0 },
              { -1, 0, 0 },
              { -1, -1, 0 }};
int N16[13][3]= {{ 0, 1, 0 },
              {-1, 2, 0 },
              {-1, 1, 0 },
              {-2, 1, 0 },
              {-1, 0, 0 },
              {-2,-1, 0 },
              {-1,-1, 0 },
              {-1,-2, 0 }};

/*-- 3D --*/
int N6[13][3] = {{ 1, 0, 0},
              { 0, 1, 0},
              { 0, 0, 1}};
int N18[13][3]= {{ 1, 0, 0},
                 { 1, 1, 0},
                 { 0, 1, 0},
                 {-1, 1, 0},
                 { 0, 0, 1},
                 { 1, 0, 1},
                 { 0, 1, 1},
                 {-1, 0, 1},
                 { 0,-1, 1}};
int N26[13][3]= {{ 1, 0, 0},
                 { 0, 1, 0},
                 { 1, 1, 0},
                 {-1, 1, 0}, /* 4 */
                 { 0, 0, 1},
                 { 1, 0, 1},
                 { 1, 1, 1},
                 { 0, 1, 1},
                 {-1, 1, 1}, /* 9 */
                 {-1, 0, 1},
                 {-1,-1, 1},
                 { 0,-1, 1}, /* 12 */
                 { 1,-1, 1}};


/* undirected neighbor edges */
static uint n_edges(int e[3], uint dims[3], uint **P, uint **Q)
{
    uint cnt = (dims[0]-abs(e[0]))*(dims[1]-abs(e[1]))*(dims[2]-abs(e[2]));//绝对值
    MALLOC_SET(uint, *P, cnt);
    MALLOC_SET(uint, *Q, cnt);

    /*- internal -*/
    const int R = e[0], C = e[1], D = e[2];
    int r_min = (R < 0)? -R : 0;
    int r_max = dims[0] - ((R < 0)? 0 : R);
    int c_min = (C < 0)? -C : 0;
    int c_max = dims[1] - ((C < 0)? 0 : C);
    int d_min = (D < 0)? -D : 0;
    int d_max = dims[2] - ((D < 0)? 0 : D);
    for (int d = d_min, i =0; d < d_max; d++) {
        for (int r = r_min; r < r_max; r++) {
            for (int c = c_min; c < c_max; c++, i++) {
                (*P)[i] = sub2ind(r, c, d) + 1;
                (*Q)[i] = sub2ind(r+R, c+C, d+D) + 1;
            }
        }
    }

    return cnt;
}
static bool usage(const char *msg)
{
    mexPrintf("error: %s\n"
              "\n"
              "usage: \n"
              " map = chop(dims, N, label_cnt, R, B, V)\n"
              "   map - labeling\n"
			  "   dims=[r c d]\n"
              "   labels=[r c d]\n"
              "   N from {4,8,26} in 2D or {6,18,26} in 3D\n"
              "   label_cnt - number of labels to use\n"
              "   R - regional term function (given label)\n"
              "   B - shared boundary term function\n"
              "   V - relative label-to-label penalties\n",
              msg);
    return false;
}

/* 1 iff scalar */
static int is_scalar(const mxArray *m)
{
    const int dims_cnt = mxGetNumberOfDimensions(m);
    const mwSize *dims = mxGetDimensions(m);
    return dims_cnt == 2 && dims[0] == 1 && dims[1] == 1;
}


static bool check_usage(int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[])
{

    /* number of parameters */
	if (nrhs != 7)
        return usage("wrong number of inputs");
    /* number of returns */
    if (nlhs != 1)
        return usage("wrong number of outputs");


    /* dims */
    int dims_cnt = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dims = mxGetDimensions(prhs[0]);
    if (dims_cnt != 2 || !(dims[0] == 1 && (dims[1] == 2 || dims[1] == 3)))
        return usage("dims must be row vector: [r c] or [r c d]");

    const double *dim = mxGetPr(prhs[0]);
    int is_3D = (dims[1] == 3) && dim[2] > 1;
    if (!(dim[0] > 1 && dim[1] > 1 && (!is_3D || dim[2] > 1)))
        return usage("dim [r c] or [r c d] for r,c at least 2, d positive");




    /* N */
    if (!is_scalar(prhs[1]))
        return usage("N must be a scalar");
    uint N = (uint)mxGetScalar(prhs[1]);
    if (is_3D && !(N == 6 || N == 18 || N == 26))
        return usage("only neighborhoods of 6,18,26 supported in 3D");
    else if (!is_3D && !(N == 4 || N == 8 || N == 16))
        return usage("only neighborhoods of 4,8,16 supported in 2D");


    /* label_cnt */
    if (!is_scalar(prhs[2]))
        return usage("label_cnt must be scalar");
    if (mxGetScalar(prhs[2]) <= 1)
        return usage("must have at least two labels");

    /* R */
    if (mxGetClassID(prhs[3]) != mxFUNCTION_CLASS)//函数类型
        return usage("R must be a function handle");

    /* B */
    if (mxGetClassID(prhs[4]) != mxFUNCTION_CLASS)
        return usage("B must be a function handle");
	/* H */	
    if ( mxGetClassID(prhs[5]) != mxSTRUCT_CLASS )
        mexErrMsgIdAndTxt("robustpn:inputs","hop must be a struct array");
    // expecting HOP_N_OF_FIELDS fieds
    if ( mxGetNumberOfFields(prhs[5]) != HOP_N_OF_FIELDS )
        mexErrMsgIdAndTxt("robustpn:inputs","hop must have %d fields", HOP_N_OF_FIELDS);
    

    /* V */
    if (!mxIsDouble(prhs[6]))
        return usage("V must be of type double");
    if (mxGetNumberOfDimensions(prhs[6]) != 2)
        return usage("V must be a matrix of n-by-n");
    const mwSize *V_dims = mxGetDimensions(prhs[6]);
    if (V_dims[0] != V_dims[1])
        return usage("V must be of size label_cnt x label_cnt");

    return true; /* passed inspection */
}

static mxArray *m_R_fn;
static mxArray *m_R_w;
static double *V;
static mxArray *m_label;
static double *R_fn()
{
    mxArray *prhs[] = {m_R_fn};
    mexCallMATLAB(1, &m_R_w, 1, prhs, "feval");
    ASSERT(mxIsDouble(m_R_w));
    return mxGetPr(m_R_w);
}


static mxArray *m_B_fn;
static mxArray *m_e;
static mxArray *m_B_w;
static double *B_fn(double *e, uint *P, uint *Q, uint cnt)
{
    mxArray *m_P = mxCreateNumericMatrix(1, cnt, mxUINT32_CLASS, mxREAL);
    mxArray *m_Q = mxCreateNumericMatrix(1, cnt, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(m_e), e, 3*sizeof(*e)); /* OPT: setdata? */
    memcpy(mxGetData(m_P), P, cnt*sizeof(*P));
    memcpy(mxGetData(m_Q), Q, cnt*sizeof(*Q));
    mxArray *prhs[] = { m_B_fn, m_e, m_P, m_Q };
    mexCallMATLAB(1, &m_B_w, 4, prhs, "feval");
    ASSERT(mxIsDouble(m_B_w));
    mxDestroyArray(m_P);
    mxDestroyArray(m_Q);
    return mxGetPr(m_B_w);
}


/*下面这个mexFunction的目的是使MATLAB知道如何调用函数*/ 
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{    /* nlhs是MATLAB命令行方式下输出参数的个数； 
     *plhs[]是MATLAB命令行方式下的输出参数； 
     nrhs是MATLAB命令行方式下输入参数的个数； 
     *prhs[]是MATLAB命令行方式下的输入参数； */ 

    mexSetTrapFlag(1); /* bail to Matlab prompt on mexCallMATLAB error */

    /*-- check command arguments --*/
    if (check_usage(nlhs, plhs, nrhs, prhs) == 0)
        return;

    /*-- input --*/
    const mwSize *dim = mxGetDimensions(prhs[0]);//该函数返回array_ptr各维的元素个数保存在一个int数组中返回
    int is_3D = (dim[1] == 3);
    double *d_dims = mxGetPr(prhs[0]);
    uint dims[3] = { d_dims[0], d_dims[1], (is_3D? d_dims[2] : 1) };
    mwSize dims_m[3] = { d_dims[0], d_dims[1], (is_3D? d_dims[2] : 1) };

    uint N = mxGetScalar(prhs[1]);
    uint label_cnt = mxGetScalar(prhs[2]);
	/*--所有 MATLAB的 C API的数据类型是 mxArray ， MATLAB 提供了一组 C  API 来进行数据操作。
    mxArray 是一种包含多种类型的数据，可以是数值，字符， cell 或者是结构体；数据类型可以是标量，矩阵,STRUCT或者是 Array 。 --*/
    m_R_fn = (mxArray *)prhs[3];
    m_B_fn = (mxArray *)prhs[4];
    V = mxGetPr(prhs[6]);


    /*-- weight calculation call backs --*/
    //m_label = mxCreateDoubleScalar(-1);//用n来初始化生成的数组,其实标量就是一个1*1的数组
    //label = mxGetPr(m_label);//获得指向label的指针
    //MALLOC_SET(mxArray *, m_R_w, label_cnt);//dest = (ty *)JM_MALLOC( (cnt) * sizeof(dest[0]) );


    /* compute neighbor weights */
    m_e = mxCreateNumericMatrix(3, 1, mxDOUBLE_CLASS, mxREAL);

	int nVar, nPair, nHigher, ii;
    int hop_fields_indices[HOP_N_OF_FIELDS];
	nVar = dims[0]*dims[1]*dims[2];
	nPair = 4*dims[0]*dims[1]-3*dims[0]-3*dims[1]+2;
	nHigher = mxGetNumberOfElements(prhs[5]);

	// chack that we have the right fields
   for (ii = 0; ii < HOP_N_OF_FIELDS ; ii++ ) {
        hop_fields_indices[ii] = mxGetFieldNumber(prhs[5], HOP_FIELDS[ii]);
        if ( hop_fields_indices[ii] < 0 )
           mexErrMsgIdAndTxt("robustpn:inputs","hop is missing %s field", HOP_FIELDS[ii]);
    }
	Energy<double> *energy = new Energy<double>(label_cnt, nVar, nPair,nHigher);
    
    // Add the HO potentials
/*      hop - higher order potential array of structs with (#higher) entries, each entry:
 *      .ind - indices of nodes belonging to this hop
 *      .w - weights w_i for each participating node
 *      .gamma - #labels + 1 entries for gamma_1..gamma_max
 *      .Q - truncation value for this potential (assumes one Q for all labels)
 */
    mxArray *xind, *xw, *xgamma, *xQ;
    int * ind, n;
    double* w;
    double* gamma;
    double Q;
    for ( ii = 0 ; ii < nHigher; ii++ ) {
        xind = mxGetFieldByNumber(prhs[5], ii, hop_fields_indices[0]);
        n = mxGetNumberOfElements(xind);
        ind = new int[n]; // allocation for energy
        GetArr(xind, ind, -1); // bias = -1 convert from 1-ind of matlab to 0-ind of C
        
        xw = mxGetFieldByNumber(prhs[5], ii, hop_fields_indices[1]);
        if ( mxGetNumberOfElements(xw) != n ) {
            delete energy;
            delete[] ind;
            mexErrMsgIdAndTxt("robustpn:inputs","hop %d: number of indices is different than number of weights", ii);
        }
        w = new double[n]; // allocation for energy
        GetArr(xw, w);
        
        xgamma = mxGetFieldByNumber(prhs[5], ii, hop_fields_indices[2]);
        if ( mxGetNumberOfElements(xgamma) != label_cnt+1 ) {
            delete energy;
            delete[] ind;
            delete[] w;
            mexErrMsgIdAndTxt("robustpn:inputs","hop %d: must have exactly %d gamma values", ii, label_cnt+1);
        }
        gamma = new double[label_cnt+1];
        GetArr(xgamma, gamma);
        
        xQ = mxGetFieldByNumber(prhs[5], ii, hop_fields_indices[3]);
        Q = (double)mxGetScalar(xQ);
        if (  energy->SetOneHOP(n, ind, w, gamma, Q)< 0 ) {

			mexErrMsgIdAndTxt("robustpn:inputs","failed to load hop %d",ii);
			delete energy;
            delete[] gamma;
        }
        delete[] gamma; // this array is being allocated inside energy
        // mexPrintf("Done reading hop(%d) / %d\n", ii, nHigher);
	}

    /*-- 调用graph cut --*/
    uint *map = graphchop(dims, N, label_cnt, R_fn, B_fn, V,energy);

	AExpand<double> *expand = new AExpand<double>(energy, MAX_ITER);

	plhs[0] = mxCreateNumericMatrix(1, nVar, mxINT32_CLASS, mxREAL);
    int *solution = (int*)mxGetData(plhs[0]);
	memset(solution, 0, nVar*sizeof(int));
    double ee[3];
    double E(0);
	mexPrintf("start\n");
    E = expand->minimize(solution, ee);
	mexPrintf("end\n");
	memcpy(mxGetData(plhs[0]), solution, dims[0]*dims[1]*dims[2]*sizeof(int));
    if (nlhs>1) {
        plhs[1] = mxCreateNumericMatrix(1, 4, mxGetClassID(prhs[1]), mxREAL);
        double *pE = (double*)mxGetData(plhs[1]);
        pE[0] = ee[0]; // unary energy
        pE[1] = ee[1]; // pair-wise energy
        pE[2] = ee[2]; // higher order energy
        pE[3] = E;
    }
    // de-allocate
    delete expand;
    delete energy;


    /*-- output --*/
	mexPrintf("start\n");
  // plhs[0] = mxCreateNumericArray(3, dims_m, mxUINT32_CLASS, mxREAL);
  // memcpy(mxGetData(plhs[0]), map, dims[0]*dims[1]*dims[2]*sizeof(*map));
   mexPrintf("end\n");

    /*-- clean up --*/
    mxDestroyArray(m_R_w);
    mxDestroyArray(m_B_w);
    mxDestroyArray(m_label);
    mxDestroyArray(m_e);
    JM_FREE(map);
}


// actual function according to desired termType class
template<class T>
void GetArr(const mxArray* x, T* arr, T bias)
{
    int ii, n = mxGetNumberOfElements(x);
    void *p = mxGetData(x);
    char* cp;
    unsigned char* ucp;
    short* sp;
    unsigned short* usp;
    int* ip;
    unsigned int* uip;
    int64_T *i64p;
    uint64_T *ui64p;
    double* dp;
    float* fp;
    
    switch (mxGetClassID(x)) {
        case mxCHAR_CLASS:
        case mxINT8_CLASS:    
            cp = (char*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = cp[ii] + bias;
            return;
        case mxDOUBLE_CLASS:
            dp = (double*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = dp[ii]+ bias;
            return;
        case mxSINGLE_CLASS:
            fp = (float*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = fp[ii]+ bias;
            return;
        case mxUINT8_CLASS:
            ucp = (unsigned char*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = ucp[ii]+ bias;
            return;
        case mxINT16_CLASS:
            sp = (short*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = sp[ii]+ bias;
            return;
        case mxUINT16_CLASS:
            usp = (unsigned short*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = usp[ii]+ bias;
            return;
        case mxINT32_CLASS:
            ip = (int*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = ip[ii]+ bias;
            return;
        case mxUINT32_CLASS:
            uip = (unsigned int*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = uip[ii]+ bias;
            return;
        case mxINT64_CLASS:
            i64p = (int64_T*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = i64p[ii]+ bias;
            return;
        case mxUINT64_CLASS:
            ui64p = (uint64_T*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = ui64p[ii]+ bias;
            return;
        default:
            mexErrMsgIdAndTxt("robustpn:GetArr","unsupported data type");
    }
}

template<typename termType>
uint *graphchop(uint dims[3], uint N, uint label_cnt,
                t_weight_fn R_fn, n_weight_fn B_fn, double *V,Energy<termType> *e)
{	
	
	/*-- build graph --*/
    double *R_w;
	uint node_cnt = dims[0]*dims[1]*dims[2];
	mexPrintf("graph\n");
	graph(dims ,N, label_cnt, R_fn, B_fn, V,&R_w,e);//dims[0]=128 dims[1]=192 dims[2]=1
   
    MALLOC_INIT(uint, map, node_cnt);
	mexPrintf("endgrapht\n");
    for (uint i = 0; i < node_cnt; i++)
		map[i] = 1;
    return map;
}

template<typename termType>
static void graph(uint dims[3],  uint N, uint label_cnt,
                           t_weight_fn R_fn, n_weight_fn B_fn, double *V,double **R_w,Energy<termType> *energy)
{	
    uint node_cnt = dims[0]*dims[1]*dims[2];
    /*-- terminal nodes --*/
    MALLOC_SET(double, *R_w, node_cnt*label_cnt);//R_w存储每个点的标号（0和1）
    *R_w = R_fn();
    /*-- neighboring system */
    int fam[13][3];  /* HACK */
    if (dims[2] > 1) {
        /* 3D */
        switch (N) {//相邻节点个数
        case 6:  memcpy(fam, N6, sizeof(fam)); break;
        case 18: memcpy(fam, N18, sizeof(fam)); break;
        case 26: memcpy(fam, N26, sizeof(fam)); break;
        default: assert(0 && "invalid 3D neighborhood");
        }
    } else {
        /* 2D */
        switch (N) {
        case 4:  memcpy(fam, N4, sizeof(fam)); break;
        case 8:  memcpy(fam, N8, sizeof(fam)); break;
        case 16: memcpy(fam, N16, sizeof(fam)); break;
        default: assert(0 && "invalid 2D neighborhood");
        }
    }
	int	max_npair = 0;//得到对数
    for (uint k = 0; k < N/2; k++) {
        uint *n_P, *n_Q;
        uint n_cnt = n_edges(fam[k], dims, &n_P, &n_Q);
		max_npair = max_npair + n_cnt;
        double e[3];
        for (int i = 0; i < 3; i++)
            e[i] = fam[k][i];
        B_fn(e, n_P, n_Q,n_cnt);//mex_main中作出了详细定义 weight.h中进行声明，用mex调用matlab程序
        JM_FREE(n_P); JM_FREE(n_Q);
    }
	mexPrintf("npair：%d\n",max_npair);
	int * pairs = new int[2 * max_npair]; // will be de-alocate on ~Energy
    termType* sc = new termType[max_npair]; // will be de-alocate on ~Energy
	int off = 0;
    for (uint k = 0; k < N/2; k++) {
        uint *n_P, *n_Q;
        uint n_cnt = n_edges(fam[k], dims, &n_P, &n_Q);
        double e[3];
        for (int i = 0; i < 3; i++)
          e[i] = fam[k][i];
        e[2]=16;
		/* set edge weights */
		// char enger[100];
        //sprintf(enger,"%d",n_w[2]);
        //MessageBox(NULL,enger," n_w[2] ",MB_OK);
        double *n_w = B_fn(e, n_P, n_Q,n_cnt);//mex_main中作出了详细定义 weight.h中进行声明，用mex调用matlab程序

        for (int i = 0; i < n_cnt; i++) {
            uint self = n_P[i] - 1;
            uint other = n_Q[i] - 1;
           // g->setNeighbors(self, other, n_w[i]);
			pairs[(i+off)*2] = self;
			pairs[(i+off)*2+1] = other;
			sc[i+off] = n_w[i]; 
        }
		off = off + n_cnt;
        JM_FREE(n_P); JM_FREE(n_Q);
    }
	mexPrintf("startset:%d\n",off);
    energy->SetUnaryCost( (termType*)(*R_w));
    energy->SetPairCost(pairs, sc);
    delete[] pairs; // were copied into energy
    delete[] sc;
}