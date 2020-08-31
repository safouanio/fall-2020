# include("C:/Users/safou/R projects/00-JuliaTips/SoluSL.jl")


using JLD2
using Random
using LinearAlgebra
using Statistics
using CSV
using DataFrames
using FreqTables


function q1()
#################### Q1 ###################
######## a
Random.seed!(1234)
A = 15 .*rand(10,7) .- 5;
B = 15 .*randn(10,7) .- 2;
C = hcat(A[ 1:5 , 1:5 ] , B[ 1:5 , 6:7 ])

latD = copy(vec(A))

for i in 1:length(latD)
    if latD[i] > 0
        latD[i] = 0
    end
end

D = reshape(latD,10,7)

######## b
length(A)
######## c
length(unique(D))
######## d
E = reshape(B,length(B),1)
Ev = vec(B)
######## e
F = cat(A,B,dims=3)
######## f
F = permutedims(F, [3,1,2])
######## g
G = kron(B,C)
# Julia rise dimensions issues
######## h
@save "C:/Users/safou/R projects/00-JuliaTips/matrixpractice.jld" A B C D E F G
######## i
@save "C:/Users/safou/R projects/00-JuliaTips/firstmatrix.jld" A B C D
######## j
CSV.write("C:/Users/safou/R projects/00-JuliaTips/Cmatrix.csv",DataFrame(C))
######## k
CSV.write("C:/Users/safou/R projects/00-JuliaTips/Dmatrix.dat",DataFrame(D), delim="\t")
######## l
return(A,B,C,D)
end

A,B,C,D = q1()




function q2(A,B,C)

#################### Q2 ###################
######## a

AB = Array{Float64}(undef,10,7)

for i in 1:size(A,1), j = 1:size(B,2)
    AB[i,j] = A[i,j]*B[i,j]
end

AB2 = A .*B
######## b

MaxElementC = size(C,1)*size(C,2)
ll = 0
CprimeI = Array{Float64}(undef,MaxElementC,1)
for i in 1:size(C,1), j = 1:size(C,2)
    if  -5 <= C[i,j] <= 5
        ll = ll + 1
        CprimeI[ll] = C[i,j]
    end
end

Cprime = CprimeI[1:ll]
Cprime2 = C[abs.(C) .<= 5]
######## c

N, K, T = 15169, 6, 5
X = Array{Float64}(undef,N,K,T)
for t in 1:T
    X[:,1,t]  = ones(N)
    X[:,2,t]  = [Int64( rand() >= (0.75 *(6 - t)/5) ) for i in 1:N ]
    X[:,3,t]  = (5*(t-1)) .*randn(N) .+ (15 + t + 1)
    X[:,4,t]  = (1/exp(1)) .*randn(N) .+ pi*(6 - t)/3
    X[:,5,t]  = trunc.((2.19^3) .*rand(N) .+ 12)
    X[:,6,t]  = trunc.((20*0.5^2) .*rand(N) .+ (20*0.5))
end

######## d
Beta = Array{Float64}(undef,K,T)
[Beta[:,t] = [0.25*t + 0.75, log(t), -sqrt(t), exp(t) - exp(t + 1), t, t/3] for t in 1:T]

######## e
Y = Array{Float64}(undef,N,T)
[Y[:,t] =  X[:,:,t]*Beta[:,t] + 0.36*randn(N) for t in 1:T]

end

q2(A,B,C)


function q3()
#################### Q3 ###################
######## a
dff = CSV.read("C:/Users/safou/R projects/00-JuliaTips/nlsw88.csv"; delim = ',') # in this data missing values are not ignored

df = CSV.read("C:/Users/safou/R projects/00-JuliaTips/nlsw88.csv"; delim = ',', missingstring = "NA")
@save "C:/Users/safou/R projects/00-JuliaTips/nlsw88.jld" df

######## b
per_nm = mean(df.never_married)*100
per_cg = mean(df.collgrad)*100

######## c
per_race = prop(freqtable(df.race)).*100
######## d

sumstdf = describe(df, :mean, :median, :std, :min, :max, :nunique, :q75, :q25)
size(sumstdf)

m = size(sumstdf)[2]-2 # last col for interquartile computing
summarystats = Array{Any}(undef,size(sumstdf)[1],m)

#initiating matrix with the discribe df value and replacing nothing values by missing which allows computing

for i in 2:m
    summarystats[:,i-1] = replace(sumstdf[:,i], nothing => missing)
end

#computing interquartile
summarystats[:,m] = replace(sumstdf[:, m + 1], nothing => missing) - replace(sumstdf[:, m + 2], nothing => missing)

#for the describe(dff) command the there is 2 missing values
describe(dff)

######## e
freqtable(df.industry, df.occupation)[1:end .!=1 , 1:end .!=1 ]  # droping row and col of missing values

######## f
subdf = select(df, :wage, :occupation, :industry)
combine(groupby(subdf, All(:industry, :occupation)), :wage => mean)

end




function q4()
#################### Q4 ###################
######## a
@load "C:/Users/safou/R projects/00-JuliaTips/firstmatrix.jld"
######## b
function matrixops(A,B)
    #compute and return a tulpe element-by-element product, matrix product and sum of elements of two matrix
    if(size(A)!=size(B))
        return("inputs must have the same size")
    end



    i = A .*B
    ii = A'B
    iii = sum(A + B)
    return(i,ii,iii)
end

######## c
# in the code above

######## d
matrixops(q1()[1],q1()[2])

######## e
# in the code above

######## f
matrixops(q1()[3],q1()[4])
# it rises the error message

######## g
convert(Array,df.ttl_exp)
convert(Array,df.wage)
matrixops(convert(Array,df.ttl_exp),convert(Array,df.wage))

######## h

end
