#include <vector>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <cfloat>

using namespace std;

template <typename REAL>
vector<vector<REAL>>& build_zero_matrix(int dim) {
    auto matr = new vector<vector<REAL>>(dim);

    for(int i = 0; i < dim; i++) {
        (*matr)[i] = vector<REAL>(dim);
        for(int j = 0; j < dim; j++) {
            (*matr)[i][j] = 0;
        }
    }

    return *matr;
}

template <typename REAL>
vector<vector<REAL>>& build_unique_matrix(int dim) {
    auto matr = new vector<vector<REAL>>(dim);

    for(int i = 0; i < dim; i++) {
        (*matr)[i] = vector<REAL>(dim);
        for(int j = 0; j < dim; j++) {
            (*matr)[i][j] = 1.0/(REAL(i+j-sin(M_PI * 3 / 2 * 13)));
        }
    }

    return *matr;
}

template <typename REAL>
vector<vector<REAL>>& scan_matr(int dim) {
	ifstream inp;
	REAL cur_elem;
	char input_t;
    vector<vector<REAL>>* matr;
    string tmp, file_name;

	cout << "How to set the matrix? (0 - from a text file; 1 - from the console; 2 - Hilbert matrix)" << endl;
	cin >> input_t;
    matr = &build_zero_matrix<REAL>(dim);

	switch(input_t) {
	    case '0': {

            //Открываем файл
            cout << "Enter the file name in CURS_21" << endl;
            cin >> file_name;
            inp.open("C:\\Users\\FireFox\\Downloads\\CURS_21\\CURS_21\\" + file_name);
            if(!inp.is_open()) {
                throw -1;
            }

            getline(inp, tmp); //пропускаем строку с размерностью

            //Считываем саму матрицу. Она пиздец по тупому у него в файл записывается.
            //Там сначала первая половина столбцов, а ПОД (!) ними вторая половина

            int step_num = dim/3+1, steps_lim;
            for(int c = 0; c < step_num+1; c++) {
                steps_lim = (dim-step_num*c >= 3) ? (c+1)*3 : dim;

                for(int i = 0; i < dim; i++) { //строки
                    for (int j = c*3; j < steps_lim; j++) {
                        inp >> (*matr)[i][j];
                    }
                }
            }
            break;
        }

        case '1': {
            cout << "Enter the matrix " << dim << "x" << dim << " by rows, with spaces:" << endl;
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    cin >> (*matr)[i][j];
                }
            }
            break;
        }

	    case '2': {
            matr = &build_unique_matrix<REAL>(dim);
            break;
        }

        default:
            matr = &build_unique_matrix<REAL>(dim);
	}

	return *matr;
}

template <typename REAL>
REAL vector_norm_1(vector<REAL> v) {
    REAL max_el = v[0];

    for(REAL el:v) {
        max_el = (abs(el) > abs(max_el)) ? abs(el) : abs(max_el);
    }
    return max_el;
}

template <typename REAL>
REAL vector_norm_2(vector<REAL> v) {
    REAL res = 0;

    for(REAL el:v) {
        res += abs(el);
    }
    return res;
}

template <typename REAL>
REAL vector_norm_3(vector<REAL> v) {
    REAL res = 0;

    for(REAL el:v) {
        res += pow(el, 2);
    }

    REAL suka = sqrt(res);

    return suka;
}

template <typename REAL>
REAL matrix_norm_2(vector<vector<REAL>>& A, int dim) {
    REAL max_col_sum = 0, col_sum;
    for(int j = 0; j < dim; j++) {
        col_sum = 0;
        for(int i = 0; i < dim; i++) {
            col_sum += abs(A[i][j]);
        }
        max_col_sum = (col_sum > max_col_sum) ? col_sum : max_col_sum;
    }

    return max_col_sum;
}

template <typename REAL>
REAL calc_matr_avg(vector<vector<REAL>>& A, int dim) {
    REAL sum = 0;
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            sum += abs(A[i][j]);
        }
    }
    return sum / pow(dim, 2);
}

template <typename REAL>
void switch_lead_rows(vector<vector<REAL>>& LU, vector<vector<REAL>>& A, int dim, int j_col, unsigned int& permtn_count, vector<REAL>& b) { //switches the row with the lead element with j_col-row
    int max_id, max_el_abs;
    vector<REAL> tmp1, tmp_A;
    REAL tmp2;

    //Сначала выберем максимальный по столбцу
    max_el_abs = abs(LU[j_col][j_col]);
    max_id = j_col;
    for(int c = j_col+1; c < dim; c++) {
        if(abs(LU[c][j_col]) > max_el_abs) {
            max_el_abs = abs(LU[c][j_col]);
            max_id = c;
        }
    }
    //Меняем местами строки в LU, исходной матрице и векторе правой части
    tmp1 = LU[max_id]; tmp2 = b[max_id]; tmp_A = A[max_id];
    LU[max_id] = LU[j_col]; b[max_id] = b[j_col]; A[max_id] = A[j_col];
    LU[j_col] = tmp1; b[j_col] = tmp2; A[j_col] = tmp_A;

    permtn_count += (max_id == j_col) ? 0 : 1;
}

template <typename REAL>
bool is_diag_suspicious(REAL diag, bool permtn, REAL matr_avg, REAL null_const) {
    if(permtn) {
        return abs(diag) <= null_const;
    }
    else {
        return abs(diag) <= null_const * matr_avg;
    }
}

template <typename REAL>
vector<vector<REAL>>& get_LU(vector<vector<REAL>>& A, int dim, bool permtn, unsigned int& ind_y, unsigned int& permtn_count,
       REAL null_const, vector<REAL>& b) {

    auto* LU = new vector<vector<REAL>>(A);
    REAL diag, r, matr_avg = calc_matr_avg(A, dim);

    for(int i = 0; i < dim-1; i++) {
        if(permtn) { //если пользователь заказал перестановку столбцов
            switch_lead_rows(*LU, A, dim, i, permtn_count, b);
        }

        //Проверяем матрицу на вырожденность
        diag = (*LU)[i][i];
        if(is_diag_suspicious(diag, permtn, matr_avg, null_const)) {
            ind_y++;
        }

        for(int k = i+1; k < dim; k++) {
            //Все поддиагональные элементы делим на диагональный
            r = (*LU)[k][i]/diag;
            (*LU)[k][i] = r;
            for(int j = i+1; j < dim; j++) {
                (*LU)[k][j] = (*LU)[k][j] - r * (*LU)[i][j];
            }
        }
    }

    return *LU;
}

template <typename REAL>
vector<REAL>& solve_lz(vector<vector<REAL>> LU, int dim, vector<REAL> b) {
    auto* z = new vector<REAL>(b);
    REAL sum;

    for(int i = 1; i < dim; i++) {
        sum = 0;
        for(int j = 0; j < i; j++) {
            sum += LU[i][j] * (*z)[j];
        }

        (*z)[i] = b[i] - sum;
    }

    return *z;
}

template <typename REAL>
vector<REAL>& solve_ux(vector<vector<REAL>> LU, int dim, vector<REAL> z) {
    auto* x = new vector<REAL>(z);
    REAL sum;
    (*x)[dim-1] /= LU[dim-1][dim-1];

    for(int i = dim-2; i >= 0; i--) {
        sum = 0;
        for(int j = i+1; j < dim; j++) {
            sum += LU[i][j] * (*x)[j];
        }

        (*x)[i] = (z[i] - sum) / LU[i][i];
    }

    return *x;
}

template<typename REAL>
vector<REAL>& solve_lu_full(vector<vector<REAL>> LU, int dim, vector<REAL> b) {
    vector<REAL> z = solve_lz(LU, dim, b);
    auto* x = new vector<REAL>(dim);

    *x = solve_ux(LU, dim, z);
    return *x;
}

template <typename REAL>
REAL scalar_multiply(const vector<REAL> &v, const vector<REAL> &u, int dim) {
    REAL res = 0;
    for(int i = 0; i < dim; i++) {
        res += v[i] * u[i];
    }
    return res;
}

template <typename REAL>
vector<REAL> vector_mult_by_scalar(const vector<REAL> &v, REAL scalar, int dim) {
    vector<REAL> u(dim);
    for(int i = 0; i < dim; i++) {
        u[i] = v[i] * scalar;
    }
    return u;
}

template <typename REAL>
vector<REAL> sum_vectors(const vector<REAL> &v, const vector<REAL> &u, int dim) {
    vector<REAL> res(dim);
    for(int i = 0; i < dim; i++) {
        res[i] = v[i] + u[i];
    }
    return res;
}

template <typename REAL>
vector<REAL>& matr_mult_by_vector(vector<vector<REAL>>& A, vector<REAL>& x, int dim) {
    auto* res = new vector<REAL>(dim);
    for(int i = 0; i < dim; i++) {
        (*res)[i] = scalar_multiply(A[i], x, dim);
    }

    return *res;
}

template <typename REAL>
REAL calc_residual(vector<vector<REAL>>& A, int dim, vector<REAL>& x, vector<REAL>& b, int vec_norm_type) {
    auto* resdl = new vector<REAL>(dim);
    vector<REAL> b_ = vector_mult_by_scalar<REAL>(b, -1.0, dim);
    vector<REAL> a_ = matr_mult_by_vector<REAL>(A, x, dim);
    REAL tmp1, tmp2;

    *resdl = sum_vectors<REAL>(a_, b_, dim);

    switch (vec_norm_type) {
        case 1:
            tmp1 = vector_norm_1(*resdl),
            tmp2 = vector_norm_1(b_);
            break;
        case 2:
            tmp1 = vector_norm_2(*resdl),
            tmp2 = vector_norm_2(b_);
            break;
        case 3:
            tmp1 = vector_norm_3(*resdl),
            tmp2 = vector_norm_3(b_);
            break;
        default:
            tmp1 = vector_norm_3(*resdl),
            tmp2 = vector_norm_3(b_);
    }


    return tmp1 / tmp2;
}

template <typename REAL>
vector<REAL>& calc_vec_residual(vector<vector<REAL>>& A, int dim, vector<REAL>& x, vector<REAL>& b) {
    auto* resdl = new vector<REAL>(dim);
    //b-Ax
    vector<REAL> a = matr_mult_by_vector<REAL>(A, x, dim),
                 a_ = vector_mult_by_scalar<REAL>(a, -1.0, dim);
                 *resdl = sum_vectors<REAL>(a_, b, dim);

    return *resdl;
}

template <typename REAL>
REAL calc_rel_resld(vector<REAL>& x, int dim, int vec_norm_type) {
    vector<REAL> reltv_soltn(dim), reltv_substr(dim);
    REAL reltv_tmp1, reltv_tmp2, reltv_resdl;

    for(int i = 0; i < dim; i++) { reltv_soltn[i] = -(i + 1); }
    reltv_substr = sum_vectors(x, reltv_soltn, dim); //x-y
    switch (vec_norm_type) {
        case 1:
            reltv_tmp1 = vector_norm_1(reltv_substr);
            reltv_tmp2 = vector_norm_1(reltv_soltn);
            break;
        case 2:
            reltv_tmp1 = vector_norm_2(reltv_substr);
            reltv_tmp2 = vector_norm_2(reltv_soltn);
            break;
        case 3:
            reltv_tmp1 = vector_norm_3(reltv_substr);
            reltv_tmp2 = vector_norm_3(reltv_soltn);
            break;
        default:
            reltv_tmp1 = vector_norm_3(reltv_substr);
            reltv_tmp2 = vector_norm_3(reltv_soltn);
    }

    reltv_resdl = reltv_tmp1 / reltv_tmp2;

    return reltv_resdl;
}

template <typename REAL>
vector<vector<REAL>>& get_L(vector<vector<REAL>>& LU, int dim) {
    auto* L = new vector<vector<REAL>>(LU);
    for(int i = 0; i < dim-1; i++) {
        (*L)[i][i] = 1.0;
        for(int j = i+1; j < dim; j++) {
            (*L)[i][j] = 0;
        }
    }
    (*L)[dim-1][dim-1] = 1.0;

    return *L;
}

template <typename REAL>
vector<vector<REAL>>& get_U(vector<vector<REAL>>& LU, int dim) {
    auto* U = new vector<vector<REAL>>(LU);
    for(int i = 1; i < dim; i++) {
        for(int j = i-1; j >= 0; j--) {
            (*U)[i][j] = 0;
        }
    }

    return *U;
}

