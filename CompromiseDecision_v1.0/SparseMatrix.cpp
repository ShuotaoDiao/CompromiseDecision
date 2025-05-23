//
//  SparseMatrix.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/17/25.
//

#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix() { // default constructor
    
}

// constructor input: std::vector<std::vector<double>>
SparseMatrix::SparseMatrix(const std::vector<std::vector<double>>& denseMatrix) {
    row_length = denseMatrix.size();
    col_length = denseMatrix[0].size();
    cbeg.push_back(0);
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        clen.push_back(0);
        for (int row_idx = 0; row_idx < row_length; ++row_idx) {
            if (denseMatrix[row_idx][col_idx] > precision || denseMatrix[row_idx][col_idx] < -precision) {
                clen[col_idx] += 1;
                rind.push_back(row_idx);
                val.push_back(denseMatrix[row_idx][col_idx]);
            }
        }
        cbeg.push_back(cbeg[col_idx] + clen[col_idx]);
    }
}

// copy constructor
SparseMatrix::SparseMatrix(const SparseMatrix& sm) {
    row_length = sm.row_length;
    col_length = sm.col_length;
    cbeg = sm.cbeg;
    clen = sm.clen;
    rind = sm.rind;
    val = sm.val;
}

void SparseMatrix::operator=(const SparseMatrix& sm) {
    row_length = sm.row_length;
    col_length = sm.col_length;
    cbeg = sm.cbeg;
    clen = sm.clen;
    rind = sm.rind;
    val = sm.val;
}

// reset matrix based on the input matrix
void SparseMatrix::setSparseMatrix(const std::vector<std::vector<double>>& denseMatrix) {
    row_length = denseMatrix.size();
    col_length = denseMatrix[0].size();
    cbeg.push_back(0);
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        clen.push_back(0);
        for (int row_idx = 0; row_idx < row_length; ++row_idx) {
            if (denseMatrix[row_idx][col_idx] > precision || denseMatrix[row_idx][col_idx] < -precision) {
                clen[col_idx] += 1;
                rind.push_back(row_idx);
                val.push_back(denseMatrix[row_idx][col_idx]);
            }
        }
        cbeg.push_back(cbeg[col_idx] + clen[col_idx]);
    }
}

// fast right multiply, assume that w is a zero vector; x' * M
SparseVector SparseMatrix::fast_rightMultiply(SparseVector& sv, std::vector<double>& w) {
    // scatter
    // vector
    sv.scatter(w);
    // matrix multiplication
    SparseVector sv_res;
    sv_res.setLen(col_length);
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        int beg_idx = cbeg[col_idx]; // obtain the beginning index for each column
        double rval = 0; // initialize row value for the object
        for (int idx = beg_idx; idx < clen[col_idx] + beg_idx; ++idx) {
            if (w[rind[idx]] > precision || w[rind[idx]] < -precision) {
                rval += w[rind[idx]] * val[idx];
            }
        }
        if (rval > precision || rval < -precision) { // rval is not 0
            sv_res.insert(col_idx, rval);
        }
    }
    // erase w
    sv.eraseWorkVector(w);
    return sv_res;
}

// sparse matrix right times dense vector, return dense vector
std::vector<double> SparseMatrix::fast_rightMultiply(std::vector<double> v) {
    std::vector<double> res;
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        int beg_idx = cbeg[col_idx]; // obtain the beginning index for each column
        double rval = 0; // initialize row value for the object
        for (int idx = beg_idx; idx < clen[col_idx] + beg_idx; ++idx) {
            rval += v[rind[idx]] * val[idx];
        }
        res.push_back(rval);
    }
    return res;
}

// fast left multiply, assume that v and w are zero vectors; M * x
SparseVector SparseMatrix::fast_leftMultiply(SparseVector& sv, std::vector<double>& w, std::vector<double>& v) {
    // scatter
    // vector
    sv.scatter(w);
    // matrix multiplication
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        if (w[col_idx] > precision || w[col_idx] < -precision) {
            int beg_idx = cbeg[col_idx]; // obtain the beginning index for each column
            for (int idx = beg_idx; idx < clen[col_idx] + beg_idx; ++idx) {
                v[rind[idx]] += w[col_idx] * val[idx];
            }
        }
    }
    // erase w
    sv.eraseWorkVector(w);
    SparseVector sv_res(v);
    // erase v
    for (int idx = 0; idx < v.size(); ++idx) {
        v[idx] = 0;
    }
    return sv_res;
}

// fast left multiply a dense vector M * v
std::vector<double> SparseMatrix::fast_leftMultiply(std::vector<double> v) {
    std::vector<double> res(row_length, 0.0);
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        int beg_idx = cbeg[col_idx]; // obtain the beginning index for each column
        for (int idx = beg_idx; idx < clen[col_idx] + beg_idx; ++idx) {
            res[rind[idx]] += v[col_idx] * val[idx];
        }
    }
    return res;
}

// insert
void SparseMatrix::insert(int row_idx, int col_idx, double value) {
    if (cur_col_idx < col_idx) {
        // new column
        if (cur_col_idx == -1) { // first column
            // update beginning index
            cbeg.push_back(0);
        }
        else {
            // update beginning index
            cbeg.push_back(cbeg[cur_col_idx] + clen[cur_col_idx]);
        }
        clen.push_back(0);
        cur_col_idx += 1;
        while (cur_col_idx < col_idx) {
            clen.push_back(0);
            cbeg.push_back(cbeg[cur_col_idx] + clen[cur_col_idx]);
            cur_col_idx += 1;
        }
        if (value > precision || value < -precision) {
            clen[col_idx] += 1;
            rind.push_back(row_idx);
            val.push_back(value);
        }
    }
    else if (cur_col_idx == col_idx) {
        if (value > precision || value < -precision) {
            clen[col_idx] += 1;
            rind.push_back(row_idx);
            val.push_back(value);
        }
    }
    else {
        throw std::out_of_range("SparseMatrix::insert(int row_idx, int col_idx, double value) Error: cur_col_idx > col_idx.\n");
    }
}

void SparseMatrix::insert_end(int input_row_len, int input_col_len) {
    row_length = input_row_len;
    col_length = input_col_len;
    cbeg.push_back(cbeg[cur_col_idx] + clen[cur_col_idx]);
}

// convert to dense matrix
std::vector<std::vector<double>> SparseMatrix::scatter() {
    std::vector<std::vector<double>> res;
    for (int row_idx = 0; row_idx < row_length; ++row_idx) {
        std::vector<double> cur_row;
        for (int col_idx = 0; col_idx < col_length; ++col_idx) {
            cur_row.push_back(0);
        }
        res.push_back(cur_row);
    }
    for (int col_idx = 0; col_idx < col_length; ++col_idx) {
        int beg_idx = cbeg[col_idx]; // obtain the beginning index for each column
        for (int idx = beg_idx; idx < clen[col_idx] + beg_idx; ++idx) {
            res[rind[idx]][col_idx] = val[idx];
        }
    }
    return res;
}
// *******************************************
// Functions for traversal
// get row
int SparseMatrix::getRow(int idx) const{
    if (idx >= rind.size()) {
        throw std::out_of_range("SparseMatrix::getRow(int idx) Error: Index is out of range.\n");
    }
    return rind[idx];
}

// get entry value
double SparseMatrix::getVal(int idx) const{
    if (idx >= val.size()) {
        throw std::out_of_range("SparseMatrix::getVal(int idx) Error: Index is out of range.\n");
    }
    return val[idx];
}

// get number of nonzero entries in each column
int SparseMatrix::getClen(int idx) const {
    if (idx > clen.size()) {
        throw std::out_of_range("SparseMatrix::getClen(int idx) Error: Index is out of range.\n");
    }
    return clen[idx];
}

// get the beginning index
int SparseMatrix::getCbeg(int idx) const {
    if (idx >= cbeg.size()) {
        throw std::out_of_range("SparseMatrix::getCbeg(int idx) Error: Index is out of range.\n");
    }
    return cbeg[idx];
}

// get row length
long SparseMatrix::getRowLength() const {
    return row_length;
}

// get column length
long SparseMatrix::getColLength() const {
    return col_length;
}

// get number of nonzero entries
long SparseMatrix::getNonZeroEntries() const {
    return val.size();
}
// *******************************************
