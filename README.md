# Better Size Estimation for Sparse Matrix Products

This project implements a fast, probabilistic algorithm for estimating the number of non-zero entries in the product of two sparse boolean matrices. The algorithm is based on the 2011 paper [*Better size estimation for sparse matrix products*](https://arxiv.org/abs/1006.4173) by Amossen, Campagna, and Pagh.

## Overview

Given two sparse relations `R1(a, b)` and `R2(b, c)`, the goal is to estimate the size of the join-project `π_ac(R1 ⨝ R2)` — equivalent to the number of non-zero entries in the boolean matrix product. 

## Usage

Please refer to main.cpp for basic usage. 
