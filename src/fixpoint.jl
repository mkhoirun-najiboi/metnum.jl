"""
    fixpoint(g,p0)

adalah fungsi untuk mencari solusi dari `x=g(x)` dengan nilai tebakan awal `p0`

# Examples
```jl
julia> g(x) = exp(-x);

julia> pn, flag, M = fixpoint(g,0.5);

julia> pn
0.5671432633594872

julia> flag
0

julia> M
27×3 Array{Float64,2}:
  0.0  0.5       NaN
  1.0  0.606531    0.106531
  2.0  0.545239    0.0612914
  3.0  0.579703    0.0344639
  ⋮
 23.0  0.567143    4.09741e-7
 24.0  0.567143    2.32382e-7
 25.0  0.567143    1.31794e-7
 26.0  0.567143    7.4746e-8
```
return solusi `pn` dengan `flag` bernilai 0 jika metode bisection berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, pk, f(pk)]`

"""
function fixpoint(g,p0)
    # Definisikan nilai toleransi, maksimum iterasi dan tebakan awal
    delta = 10^-7;
    maxi = 100;
    flag = 1;
    pn = p0
    M = [0 pn NaN];
    for n = 2:maxi
        pn1 = pn;
        pn = g(pn1);
        err = abs(pn-pn1);
        relerr = err/(abs(pn)+eps());
        M = [M; [n-1 pn err]]
        if (err<delta) || (relerr<delta)
            flag = 0; break
        end
    end
    return pn, flag, M
end
