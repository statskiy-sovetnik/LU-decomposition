#include <iostream>
#include "gauss_lib.h"
#include <clocale>
#include <windows.h>
#include <iomanip>

void print_intro() {
    cout << "Gladkikh Ivan Evgenievich. MMF NSU 2018-2022. group 18121\n" << "__________________________________\n" << endl;
}

template <typename REAL>
void set_type_null_const(REAL& null_const);

template<> void set_type_null_const<float>(float &null_const) {
    null_const = FLT_MIN;
}

template<> void set_type_null_const<double> (double &null_const) {
    null_const = DBL_MIN;
}

template <typename REAL>
vector<REAL>& scan_vector(int dim) {
    auto* b = new vector<REAL>(dim);
    for(int i = 0; i < dim; i++) {
        cin >> (*b)[i];
    }
    return *b;
}

template <typename REAL>
void print_matr(vector<vector<REAL>>& A, int dim) {
    cout << "[";
    for(vector<REAL>& v:A) {
        cout << ' ';
        for(REAL el:v) {
            cout << el << ' ';
        }
        cout << endl;
    }
    cout << "]" << endl;
}

template <typename REAL>
vector<vector<REAL>>& mult_matrices(vector<vector<REAL>>& A, vector<vector<REAL>>& B, int dim) {
    vector<vector<REAL>>* RES = &build_zero_matrix<REAL>(dim);

    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            for(int k = 0; k < dim; k++)
                (*RES)[i][j] += A[i][k] * B[k][j];
        }
    }

    return *RES;
}

template <typename REAL>
int sg(REAL value) {
    return (value >= 0) ? 1 : -1;
}

template <typename REAL>
void transpose(vector<vector<REAL>>& A, int dim) {
    vector<vector<REAL>> A_copy = A;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            A[i][j] = A_copy[j][i];
        }
    }
}

template <typename REAL>
void print_errors(vector<vector<REAL>>& A, vector<vector<REAL>>& A_inv, vector<REAL> x, vector<REAL> b, int dim, int permtn_count) {

    //Выводим перестановки
    cout << "Permutations: " << permtn_count << endl;

    bool make_itr_ref = false, calc_dif_norms = false;
    int itr_ref_steps = 0;
    REAL resdl, reltv_resdl;
    vector<REAL> vec_resdl(dim), delta_x(dim);

    do {
        //Итерационное уточнение

        for(int c = 0; c < itr_ref_steps; c++) {
            vec_resdl = calc_vec_residual(A, dim, x, b);
            delta_x = solve_lu_full(A, dim, vec_resdl);
            x = sum_vectors(x, delta_x, dim);
        }

        //Выводим невязку
        resdl = calc_residual(A, dim, x, b, 3);
        cout << "Residual: " << setprecision(18) << resdl << endl;

        //Относительная невязка
        reltv_resdl = calc_rel_resld(x, dim, 3);
        cout << "Relative residual: " << setprecision(18) << reltv_resdl << endl;

        //Невязки в других нормах
        cout << "Residuals above were calculated with norm_3. Calculate with different norms? (0-no; 1-yes)" << endl;
        cin >> calc_dif_norms;
        if(calc_dif_norms) {
            cout << "Residual (norm_1): " << calc_residual(A, dim, x, b, 1) << endl;
            cout << "Residual (norm_2): " << calc_residual(A, dim, x, b, 2) << endl;
            cout << "Relative residual (norm_1): " << calc_rel_resld(x, dim, 1) << endl;
            cout << "Relative residual (norm_2): " << calc_rel_resld(x, dim, 2) << endl;
        }

        cout << "Make iterative refinement? (0-no; 1-yes)" << endl;
        cin >> make_itr_ref;
        if(make_itr_ref) {
            cout << "Enter the number of iterations" << endl;
            cin >> itr_ref_steps;
        }

    } while (make_itr_ref);


    //Число обусловленности

    REAL cond_A = matrix_norm_2(A, dim) * matrix_norm_2(A_inv, dim),
         cond_resdl = resdl * cond_A;
    cout << "Cond(A) = " << cond_A << endl;
    cout << "Residual * Cond(A) = " << cond_resdl << endl;
    cout << endl;

}

template <typename REAL>
void main_LU_func(int dim) {
    bool permtn, ind_decsn, matr_rednt = false;
    unsigned int ind_y = 0, permtn_count = 0;
    REAL null_const;

    //Инициализируем нулевую константу
    set_type_null_const<REAL>(null_const);

    //Считываем матрицу и вектор правой части
    vector<vector<REAL>>& A = build_zero_matrix<REAL>(dim);
    try{
        A = scan_matr<REAL>(dim);
    }
    catch (int er) {
        throw er;
    }
    cout << "Now enter the right part vector b: " << endl;
    vector<REAL>& b = scan_vector<REAL>(dim);

    //Берем LU разложение
    cout << "Enalbe partial pivoting? (0 - no; 1 - yes)" << endl;
    cin >> permtn;
    vector<vector<REAL>>& LU = get_LU<REAL>(A, dim, permtn, ind_y, permtn_count, null_const, b);
    vector<vector<REAL>>& L = get_L(LU, dim), U = get_U(LU, dim);

    //Проверим А на вырожденность (если не было очень маленьких диагональных элементов)
    if(ind_y == 0) {
        REAL log_sum = 0;
        //int sgn = pow(-1, permtn_count);

        for(int i = 0; i < dim; i++) {
            log_sum += log10(abs(LU[i][i])); //может быть очень малым, но не нулем
            //sgn *= sg(LU[i][i]);
        }
        //log_sum *= sgn;

        if(pow(10, log_sum) <= null_const) { //матрица А вырождена
            matr_rednt = true;
        }
    }

    //Выводим L и U
    cout << "Matrix L: " << endl;
    print_matr(L, dim); cout << endl;
    cout << "Matrix U: " << endl;
    print_matr(U, dim); cout << endl;

    if(ind_y > 0 || matr_rednt) {
        cout << "The given matrix is redundant. Continue solving the system? (0 - no; 1 - yes)" << endl;
        cin >> ind_decsn;
        if(!ind_decsn) { return; }
    }

    //Решаем систему
    vector<REAL>& z = solve_lz(LU, dim, b);
    vector<REAL>& x = solve_ux(LU, dim, z);

    //Находим обратную к A
    vector<vector<REAL>> A_inv = build_zero_matrix<REAL>(dim);
    vector<REAL> z_i, e_i, b_i;

    for(int i = 0; i < dim; i++) {
        b_i = vector<REAL>(dim, 0); b_i[i] = 1;  //инициализируем вектор Кронекера
        z_i = solve_lz(LU, dim, b_i);
        e_i = solve_ux(LU, dim, z_i);
        A_inv[i] = e_i;
    }
    transpose(A_inv, dim);


    //Выводим всякие погрешности:
    print_errors(A, A_inv, x, b, dim, permtn_count);
}


int main() {

    //выбираем нужную кодировку для потоков ввода и вывода
    //setlocale(LC_ALL, "Russian");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    char data_t = ' ';
    int dim = 0;

    print_intro();

    cout << "Data type - (0 - float; 1 - double)" << endl;
    cin >> data_t;
    cout << "Enter matrix dimension: " << endl;
    cin >> dim;

    try {
        switch (data_t) {
            case '0':
                main_LU_func<float>(dim);
                break;
            case '1':
                main_LU_func<double>(dim);
                break;
            default:
                main_LU_func<double>(dim);
        }
    }
    catch(int er) {
        if(er == -1) {
            cout << "No such file" << endl;
        }
    }

    return 0;
}