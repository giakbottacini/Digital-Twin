ЛЗ
џЯ
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype

Conv2D

input"T
filter"T
output"T"
Ttype:	
2"
strides	list(int)"
use_cudnn_on_gpubool(",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 "-
data_formatstringNHWC:
NHWCNCHW" 
	dilations	list(int)

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	

MaxPool

input"T
output"T"
Ttype0:
2	"
ksize	list(int)(0"
strides	list(int)(0",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 ":
data_formatstringNHWC:
NHWCNCHWNCHW_VECT_C
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
9
Softmax
logits"T
softmax"T"
Ttype:
2
N
Squeeze

input"T
output"T"	
Ttype"
squeeze_dims	list(int)
 (
О
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring 
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
-
Tanh
x"T
y"T"
Ttype:

2

VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 "serve*2.5.02v2.5.0-rc3-213-ga4dfb8d1a718ЫЭ
z
Conv_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameConv_1/kernel
s
!Conv_1/kernel/Read/ReadVariableOpReadVariableOpConv_1/kernel*"
_output_shapes
: *
dtype0
n
Conv_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameConv_1/bias
g
Conv_1/bias/Read/ReadVariableOpReadVariableOpConv_1/bias*
_output_shapes
: *
dtype0
z
Conv_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*
shared_nameConv_2/kernel
s
!Conv_2/kernel/Read/ReadVariableOpReadVariableOpConv_2/kernel*"
_output_shapes
: @*
dtype0
n
Conv_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameConv_2/bias
g
Conv_2/bias/Read/ReadVariableOpReadVariableOpConv_2/bias*
_output_shapes
:@*
dtype0
z
Conv_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:@ *
shared_nameConv_3/kernel
s
!Conv_3/kernel/Read/ReadVariableOpReadVariableOpConv_3/kernel*"
_output_shapes
:@ *
dtype0
n
Conv_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameConv_3/bias
g
Conv_3/bias/Read/ReadVariableOpReadVariableOpConv_3/bias*
_output_shapes
: *
dtype0
y
Dense_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	@*
shared_nameDense_1/kernel
r
"Dense_1/kernel/Read/ReadVariableOpReadVariableOpDense_1/kernel*
_output_shapes
:	@*
dtype0
p
Dense_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameDense_1/bias
i
 Dense_1/bias/Read/ReadVariableOpReadVariableOpDense_1/bias*
_output_shapes
:@*
dtype0
x
Dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@*
shared_nameDense_2/kernel
q
"Dense_2/kernel/Read/ReadVariableOpReadVariableOpDense_2/kernel*
_output_shapes

:@*
dtype0
p
Dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameDense_2/bias
i
 Dense_2/bias/Read/ReadVariableOpReadVariableOpDense_2/bias*
_output_shapes
:*
dtype0

Dense_3_with_Softmax/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*,
shared_nameDense_3_with_Softmax/kernel

/Dense_3_with_Softmax/kernel/Read/ReadVariableOpReadVariableOpDense_3_with_Softmax/kernel*
_output_shapes

:*
dtype0

Dense_3_with_Softmax/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:**
shared_nameDense_3_with_Softmax/bias

-Dense_3_with_Softmax/bias/Read/ReadVariableOpReadVariableOpDense_3_with_Softmax/bias*
_output_shapes
:*
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0

Adam/Conv_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/Conv_1/kernel/m

(Adam/Conv_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Conv_1/kernel/m*"
_output_shapes
: *
dtype0
|
Adam/Conv_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/Conv_1/bias/m
u
&Adam/Conv_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/Conv_1/bias/m*
_output_shapes
: *
dtype0

Adam/Conv_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*%
shared_nameAdam/Conv_2/kernel/m

(Adam/Conv_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Conv_2/kernel/m*"
_output_shapes
: @*
dtype0
|
Adam/Conv_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*#
shared_nameAdam/Conv_2/bias/m
u
&Adam/Conv_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/Conv_2/bias/m*
_output_shapes
:@*
dtype0

Adam/Conv_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:@ *%
shared_nameAdam/Conv_3/kernel/m

(Adam/Conv_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Conv_3/kernel/m*"
_output_shapes
:@ *
dtype0
|
Adam/Conv_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/Conv_3/bias/m
u
&Adam/Conv_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/Conv_3/bias/m*
_output_shapes
: *
dtype0

Adam/Dense_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	@*&
shared_nameAdam/Dense_1/kernel/m

)Adam/Dense_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Dense_1/kernel/m*
_output_shapes
:	@*
dtype0
~
Adam/Dense_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*$
shared_nameAdam/Dense_1/bias/m
w
'Adam/Dense_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/Dense_1/bias/m*
_output_shapes
:@*
dtype0

Adam/Dense_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@*&
shared_nameAdam/Dense_2/kernel/m

)Adam/Dense_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Dense_2/kernel/m*
_output_shapes

:@*
dtype0
~
Adam/Dense_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/Dense_2/bias/m
w
'Adam/Dense_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/Dense_2/bias/m*
_output_shapes
:*
dtype0
 
"Adam/Dense_3_with_Softmax/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*3
shared_name$"Adam/Dense_3_with_Softmax/kernel/m

6Adam/Dense_3_with_Softmax/kernel/m/Read/ReadVariableOpReadVariableOp"Adam/Dense_3_with_Softmax/kernel/m*
_output_shapes

:*
dtype0

 Adam/Dense_3_with_Softmax/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*1
shared_name" Adam/Dense_3_with_Softmax/bias/m

4Adam/Dense_3_with_Softmax/bias/m/Read/ReadVariableOpReadVariableOp Adam/Dense_3_with_Softmax/bias/m*
_output_shapes
:*
dtype0

Adam/Conv_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/Conv_1/kernel/v

(Adam/Conv_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Conv_1/kernel/v*"
_output_shapes
: *
dtype0
|
Adam/Conv_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/Conv_1/bias/v
u
&Adam/Conv_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/Conv_1/bias/v*
_output_shapes
: *
dtype0

Adam/Conv_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*%
shared_nameAdam/Conv_2/kernel/v

(Adam/Conv_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Conv_2/kernel/v*"
_output_shapes
: @*
dtype0
|
Adam/Conv_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*#
shared_nameAdam/Conv_2/bias/v
u
&Adam/Conv_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/Conv_2/bias/v*
_output_shapes
:@*
dtype0

Adam/Conv_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:@ *%
shared_nameAdam/Conv_3/kernel/v

(Adam/Conv_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Conv_3/kernel/v*"
_output_shapes
:@ *
dtype0
|
Adam/Conv_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/Conv_3/bias/v
u
&Adam/Conv_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/Conv_3/bias/v*
_output_shapes
: *
dtype0

Adam/Dense_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	@*&
shared_nameAdam/Dense_1/kernel/v

)Adam/Dense_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Dense_1/kernel/v*
_output_shapes
:	@*
dtype0
~
Adam/Dense_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*$
shared_nameAdam/Dense_1/bias/v
w
'Adam/Dense_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/Dense_1/bias/v*
_output_shapes
:@*
dtype0

Adam/Dense_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@*&
shared_nameAdam/Dense_2/kernel/v

)Adam/Dense_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Dense_2/kernel/v*
_output_shapes

:@*
dtype0
~
Adam/Dense_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/Dense_2/bias/v
w
'Adam/Dense_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/Dense_2/bias/v*
_output_shapes
:*
dtype0
 
"Adam/Dense_3_with_Softmax/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*3
shared_name$"Adam/Dense_3_with_Softmax/kernel/v

6Adam/Dense_3_with_Softmax/kernel/v/Read/ReadVariableOpReadVariableOp"Adam/Dense_3_with_Softmax/kernel/v*
_output_shapes

:*
dtype0

 Adam/Dense_3_with_Softmax/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*1
shared_name" Adam/Dense_3_with_Softmax/bias/v

4Adam/Dense_3_with_Softmax/bias/v/Read/ReadVariableOpReadVariableOp Adam/Dense_3_with_Softmax/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
ѕO
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*АO
valueІOBЃO BO
Ч
layer-0
layer_with_weights-0
layer-1
layer-2
layer-3
layer_with_weights-1
layer-4
layer-5
layer-6
layer_with_weights-2
layer-7
	layer-8

layer-9
layer-10
layer_with_weights-3
layer-11
layer_with_weights-4
layer-12
layer_with_weights-5
layer-13
	optimizer
	variables
regularization_losses
trainable_variables
	keras_api

signatures
 
h

kernel
bias
	variables
regularization_losses
trainable_variables
	keras_api
R
	variables
regularization_losses
trainable_variables
	keras_api
R
	variables
 regularization_losses
!trainable_variables
"	keras_api
h

#kernel
$bias
%	variables
&regularization_losses
'trainable_variables
(	keras_api
R
)	variables
*regularization_losses
+trainable_variables
,	keras_api
R
-	variables
.regularization_losses
/trainable_variables
0	keras_api
h

1kernel
2bias
3	variables
4regularization_losses
5trainable_variables
6	keras_api
R
7	variables
8regularization_losses
9trainable_variables
:	keras_api
R
;	variables
<regularization_losses
=trainable_variables
>	keras_api
R
?	variables
@regularization_losses
Atrainable_variables
B	keras_api
h

Ckernel
Dbias
E	variables
Fregularization_losses
Gtrainable_variables
H	keras_api
h

Ikernel
Jbias
K	variables
Lregularization_losses
Mtrainable_variables
N	keras_api
h

Okernel
Pbias
Q	variables
Rregularization_losses
Strainable_variables
T	keras_api

Uiter

Vbeta_1

Wbeta_2
	XdecaymЊmЋ#mЌ$m­1mЎ2mЏCmАDmБImВJmГOmДPmЕvЖvЗ#vИ$vЙ1vК2vЛCvМDvНIvОJvПOvРPvС
V
0
1
#2
$3
14
25
C6
D7
I8
J9
O10
P11
 
V
0
1
#2
$3
14
25
C6
D7
I8
J9
O10
P11
­
Ylayer_regularization_losses
	variables
Zmetrics

[layers
regularization_losses
trainable_variables
\non_trainable_variables
]layer_metrics
 
YW
VARIABLE_VALUEConv_1/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEConv_1/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1
 

0
1
­
^layer_regularization_losses
	variables
_metrics

`layers
regularization_losses
trainable_variables
anon_trainable_variables
blayer_metrics
 
 
 
­
clayer_regularization_losses
	variables
dmetrics

elayers
regularization_losses
trainable_variables
fnon_trainable_variables
glayer_metrics
 
 
 
­
hlayer_regularization_losses
	variables
imetrics

jlayers
 regularization_losses
!trainable_variables
knon_trainable_variables
llayer_metrics
YW
VARIABLE_VALUEConv_2/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEConv_2/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE

#0
$1
 

#0
$1
­
mlayer_regularization_losses
%	variables
nmetrics

olayers
&regularization_losses
'trainable_variables
pnon_trainable_variables
qlayer_metrics
 
 
 
­
rlayer_regularization_losses
)	variables
smetrics

tlayers
*regularization_losses
+trainable_variables
unon_trainable_variables
vlayer_metrics
 
 
 
­
wlayer_regularization_losses
-	variables
xmetrics

ylayers
.regularization_losses
/trainable_variables
znon_trainable_variables
{layer_metrics
YW
VARIABLE_VALUEConv_3/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEConv_3/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE

10
21
 

10
21
Ў
|layer_regularization_losses
3	variables
}metrics

~layers
4regularization_losses
5trainable_variables
non_trainable_variables
layer_metrics
 
 
 
В
 layer_regularization_losses
7	variables
metrics
layers
8regularization_losses
9trainable_variables
non_trainable_variables
layer_metrics
 
 
 
В
 layer_regularization_losses
;	variables
metrics
layers
<regularization_losses
=trainable_variables
non_trainable_variables
layer_metrics
 
 
 
В
 layer_regularization_losses
?	variables
metrics
layers
@regularization_losses
Atrainable_variables
non_trainable_variables
layer_metrics
ZX
VARIABLE_VALUEDense_1/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEDense_1/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE

C0
D1
 

C0
D1
В
 layer_regularization_losses
E	variables
metrics
layers
Fregularization_losses
Gtrainable_variables
non_trainable_variables
layer_metrics
ZX
VARIABLE_VALUEDense_2/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEDense_2/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE

I0
J1
 

I0
J1
В
 layer_regularization_losses
K	variables
metrics
layers
Lregularization_losses
Mtrainable_variables
non_trainable_variables
layer_metrics
ge
VARIABLE_VALUEDense_3_with_Softmax/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE
ca
VARIABLE_VALUEDense_3_with_Softmax/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE

O0
P1
 

O0
P1
В
 layer_regularization_losses
Q	variables
metrics
layers
Rregularization_losses
Strainable_variables
non_trainable_variables
layer_metrics
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
 

0
 1
f
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
8

Ёtotal

Ђcount
Ѓ	variables
Є	keras_api
I

Ѕtotal

Іcount
Ї
_fn_kwargs
Ј	variables
Љ	keras_api
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

Ё0
Ђ1

Ѓ	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

Ѕ0
І1

Ј	variables
|z
VARIABLE_VALUEAdam/Conv_1/kernel/mRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/Conv_1/bias/mPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/Conv_2/kernel/mRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/Conv_2/bias/mPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/Conv_3/kernel/mRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/Conv_3/bias/mPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/Dense_1/kernel/mRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/Dense_1/bias/mPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/Dense_2/kernel/mRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/Dense_2/bias/mPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE"Adam/Dense_3_with_Softmax/kernel/mRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam/Dense_3_with_Softmax/bias/mPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/Conv_1/kernel/vRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/Conv_1/bias/vPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/Conv_2/kernel/vRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/Conv_2/bias/vPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/Conv_3/kernel/vRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/Conv_3/bias/vPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/Dense_1/kernel/vRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/Dense_1/bias/vPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/Dense_2/kernel/vRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/Dense_2/bias/vPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE"Adam/Dense_3_with_Softmax/kernel/vRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE Adam/Dense_3_with_Softmax/bias/vPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

$serving_default_Convolutional_inputsPlaceholder*,
_output_shapes
:џџџџџџџџџШ*
dtype0*!
shape:џџџџџџџџџШ

StatefulPartitionedCallStatefulPartitionedCall$serving_default_Convolutional_inputsConv_1/kernelConv_1/biasConv_2/kernelConv_2/biasConv_3/kernelConv_3/biasDense_1/kernelDense_1/biasDense_2/kernelDense_2/biasDense_3_with_Softmax/kernelDense_3_with_Softmax/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8 *-
f(R&
$__inference_signature_wrapper_920102
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 

StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename!Conv_1/kernel/Read/ReadVariableOpConv_1/bias/Read/ReadVariableOp!Conv_2/kernel/Read/ReadVariableOpConv_2/bias/Read/ReadVariableOp!Conv_3/kernel/Read/ReadVariableOpConv_3/bias/Read/ReadVariableOp"Dense_1/kernel/Read/ReadVariableOp Dense_1/bias/Read/ReadVariableOp"Dense_2/kernel/Read/ReadVariableOp Dense_2/bias/Read/ReadVariableOp/Dense_3_with_Softmax/kernel/Read/ReadVariableOp-Dense_3_with_Softmax/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp(Adam/Conv_1/kernel/m/Read/ReadVariableOp&Adam/Conv_1/bias/m/Read/ReadVariableOp(Adam/Conv_2/kernel/m/Read/ReadVariableOp&Adam/Conv_2/bias/m/Read/ReadVariableOp(Adam/Conv_3/kernel/m/Read/ReadVariableOp&Adam/Conv_3/bias/m/Read/ReadVariableOp)Adam/Dense_1/kernel/m/Read/ReadVariableOp'Adam/Dense_1/bias/m/Read/ReadVariableOp)Adam/Dense_2/kernel/m/Read/ReadVariableOp'Adam/Dense_2/bias/m/Read/ReadVariableOp6Adam/Dense_3_with_Softmax/kernel/m/Read/ReadVariableOp4Adam/Dense_3_with_Softmax/bias/m/Read/ReadVariableOp(Adam/Conv_1/kernel/v/Read/ReadVariableOp&Adam/Conv_1/bias/v/Read/ReadVariableOp(Adam/Conv_2/kernel/v/Read/ReadVariableOp&Adam/Conv_2/bias/v/Read/ReadVariableOp(Adam/Conv_3/kernel/v/Read/ReadVariableOp&Adam/Conv_3/bias/v/Read/ReadVariableOp)Adam/Dense_1/kernel/v/Read/ReadVariableOp'Adam/Dense_1/bias/v/Read/ReadVariableOp)Adam/Dense_2/kernel/v/Read/ReadVariableOp'Adam/Dense_2/bias/v/Read/ReadVariableOp6Adam/Dense_3_with_Softmax/kernel/v/Read/ReadVariableOp4Adam/Dense_3_with_Softmax/bias/v/Read/ReadVariableOpConst*9
Tin2
02.	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *(
f#R!
__inference__traced_save_921139
	
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameConv_1/kernelConv_1/biasConv_2/kernelConv_2/biasConv_3/kernelConv_3/biasDense_1/kernelDense_1/biasDense_2/kernelDense_2/biasDense_3_with_Softmax/kernelDense_3_with_Softmax/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decaytotalcounttotal_1count_1Adam/Conv_1/kernel/mAdam/Conv_1/bias/mAdam/Conv_2/kernel/mAdam/Conv_2/bias/mAdam/Conv_3/kernel/mAdam/Conv_3/bias/mAdam/Dense_1/kernel/mAdam/Dense_1/bias/mAdam/Dense_2/kernel/mAdam/Dense_2/bias/m"Adam/Dense_3_with_Softmax/kernel/m Adam/Dense_3_with_Softmax/bias/mAdam/Conv_1/kernel/vAdam/Conv_1/bias/vAdam/Conv_2/kernel/vAdam/Conv_2/bias/vAdam/Conv_3/kernel/vAdam/Conv_3/bias/vAdam/Dense_1/kernel/vAdam/Dense_1/bias/vAdam/Dense_2/kernel/vAdam/Dense_2/bias/v"Adam/Dense_3_with_Softmax/kernel/v Adam/Dense_3_with_Softmax/bias/v*8
Tin1
/2-*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *+
f&R$
"__inference__traced_restore_921281Ж
н
_
C__inference_flatten_layer_call_and_return_conditional_losses_919245

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
Ь
d
E__inference_dropout_2_layer_call_and_return_conditional_losses_920709

inputs
identityc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/ShapeИ
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџ *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2
dropout/GreaterEqual/yТ
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџ 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
њ!
џ
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_920852

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identityЂBiasAdd/ReadVariableOpЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAdda
SoftmaxSoftmaxBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Softmaxн
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulж
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mul
IdentityIdentitySoftmax:softmax:0^BiasAdd/ReadVariableOp<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
к
c
*__inference_dropout_1_layer_call_fn_920616

inputs
identityЂStatefulPartitionedCallп
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9195232
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ&@22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
Ь
d
E__inference_dropout_1_layer_call_and_return_conditional_losses_919523

inputs
identityc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/ShapeИ
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2
dropout/GreaterEqual/yТ
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџ&@2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ&@:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
ф

И
/__inference_FCN_classifier_layer_call_fn_920131

inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identityЂStatefulPartitionedCall§
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8 *S
fNRL
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_9194072
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
ф

И
/__inference_FCN_classifier_layer_call_fn_920160

inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identityЂStatefulPartitionedCall§
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8 *S
fNRL
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_9197132
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
н[

__inference__traced_save_921139
file_prefix,
(savev2_conv_1_kernel_read_readvariableop*
&savev2_conv_1_bias_read_readvariableop,
(savev2_conv_2_kernel_read_readvariableop*
&savev2_conv_2_bias_read_readvariableop,
(savev2_conv_3_kernel_read_readvariableop*
&savev2_conv_3_bias_read_readvariableop-
)savev2_dense_1_kernel_read_readvariableop+
'savev2_dense_1_bias_read_readvariableop-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop:
6savev2_dense_3_with_softmax_kernel_read_readvariableop8
4savev2_dense_3_with_softmax_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop3
/savev2_adam_conv_1_kernel_m_read_readvariableop1
-savev2_adam_conv_1_bias_m_read_readvariableop3
/savev2_adam_conv_2_kernel_m_read_readvariableop1
-savev2_adam_conv_2_bias_m_read_readvariableop3
/savev2_adam_conv_3_kernel_m_read_readvariableop1
-savev2_adam_conv_3_bias_m_read_readvariableop4
0savev2_adam_dense_1_kernel_m_read_readvariableop2
.savev2_adam_dense_1_bias_m_read_readvariableop4
0savev2_adam_dense_2_kernel_m_read_readvariableop2
.savev2_adam_dense_2_bias_m_read_readvariableopA
=savev2_adam_dense_3_with_softmax_kernel_m_read_readvariableop?
;savev2_adam_dense_3_with_softmax_bias_m_read_readvariableop3
/savev2_adam_conv_1_kernel_v_read_readvariableop1
-savev2_adam_conv_1_bias_v_read_readvariableop3
/savev2_adam_conv_2_kernel_v_read_readvariableop1
-savev2_adam_conv_2_bias_v_read_readvariableop3
/savev2_adam_conv_3_kernel_v_read_readvariableop1
-savev2_adam_conv_3_bias_v_read_readvariableop4
0savev2_adam_dense_1_kernel_v_read_readvariableop2
.savev2_adam_dense_1_bias_v_read_readvariableop4
0savev2_adam_dense_2_kernel_v_read_readvariableop2
.savev2_adam_dense_2_bias_v_read_readvariableopA
=savev2_adam_dense_3_with_softmax_kernel_v_read_readvariableop?
;savev2_adam_dense_3_with_softmax_bias_v_read_readvariableop
savev2_const

identity_1ЂMergeV2Checkpoints
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Constl
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part2	
Const_1
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shardІ
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*
valueB-B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_namesт
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*m
valuedBb-B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slicesл
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0(savev2_conv_1_kernel_read_readvariableop&savev2_conv_1_bias_read_readvariableop(savev2_conv_2_kernel_read_readvariableop&savev2_conv_2_bias_read_readvariableop(savev2_conv_3_kernel_read_readvariableop&savev2_conv_3_bias_read_readvariableop)savev2_dense_1_kernel_read_readvariableop'savev2_dense_1_bias_read_readvariableop)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop6savev2_dense_3_with_softmax_kernel_read_readvariableop4savev2_dense_3_with_softmax_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop/savev2_adam_conv_1_kernel_m_read_readvariableop-savev2_adam_conv_1_bias_m_read_readvariableop/savev2_adam_conv_2_kernel_m_read_readvariableop-savev2_adam_conv_2_bias_m_read_readvariableop/savev2_adam_conv_3_kernel_m_read_readvariableop-savev2_adam_conv_3_bias_m_read_readvariableop0savev2_adam_dense_1_kernel_m_read_readvariableop.savev2_adam_dense_1_bias_m_read_readvariableop0savev2_adam_dense_2_kernel_m_read_readvariableop.savev2_adam_dense_2_bias_m_read_readvariableop=savev2_adam_dense_3_with_softmax_kernel_m_read_readvariableop;savev2_adam_dense_3_with_softmax_bias_m_read_readvariableop/savev2_adam_conv_1_kernel_v_read_readvariableop-savev2_adam_conv_1_bias_v_read_readvariableop/savev2_adam_conv_2_kernel_v_read_readvariableop-savev2_adam_conv_2_bias_v_read_readvariableop/savev2_adam_conv_3_kernel_v_read_readvariableop-savev2_adam_conv_3_bias_v_read_readvariableop0savev2_adam_dense_1_kernel_v_read_readvariableop.savev2_adam_dense_1_bias_v_read_readvariableop0savev2_adam_dense_2_kernel_v_read_readvariableop.savev2_adam_dense_2_bias_v_read_readvariableop=savev2_adam_dense_3_with_softmax_kernel_v_read_readvariableop;savev2_adam_dense_3_with_softmax_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *;
dtypes1
/2-	2
SaveV2К
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixesЁ
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identitym

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*№
_input_shapesо
л: : : : @:@:@ : :	@:@:@:::: : : : : : : : : : : @:@:@ : :	@:@:@:::: : : @:@:@ : :	@:@:@:::: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:($
"
_output_shapes
: : 

_output_shapes
: :($
"
_output_shapes
: @: 

_output_shapes
:@:($
"
_output_shapes
:@ : 

_output_shapes
: :%!

_output_shapes
:	@: 

_output_shapes
:@:$	 

_output_shapes

:@: 


_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :($
"
_output_shapes
: : 

_output_shapes
: :($
"
_output_shapes
: @: 

_output_shapes
:@:($
"
_output_shapes
:@ : 

_output_shapes
: :%!

_output_shapes
:	@: 

_output_shapes
:@:$ 

_output_shapes

:@: 

_output_shapes
::$ 

_output_shapes

::  

_output_shapes
::(!$
"
_output_shapes
: : "

_output_shapes
: :(#$
"
_output_shapes
: @: $

_output_shapes
:@:(%$
"
_output_shapes
:@ : &

_output_shapes
: :%'!

_output_shapes
:	@: (

_output_shapes
:@:$) 

_output_shapes

:@: *

_output_shapes
::$+ 

_output_shapes

:: ,

_output_shapes
::-

_output_shapes
: 

c
E__inference_dropout_2_layer_call_and_return_conditional_losses_920697

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
и$
ѓ
B__inference_Conv_1_layer_call_and_return_conditional_losses_920530

inputsA
+conv1d_expanddims_1_readvariableop_resource: -
biasadd_readvariableop_resource: 
identityЂBiasAdd/ReadVariableOpЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂ"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
conv1d/ExpandDims/dim
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџШ2
conv1d/ExpandDimsИ
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dimЗ
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
conv1d/ExpandDims_1И
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:џџџџџџџџџА *
paddingVALID*
strides
2
conv1d
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:џџџџџџџџџА *
squeeze_dims

§џџџџџџџџ2
conv1d/Squeeze
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:џџџџџџџџџА 2	
BiasAdd]
TanhTanhBiasAdd:output:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
Tanhв
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulК
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:џџџџџџџџџА 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:џџџџџџџџџШ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
Т
Б
__inference_loss_fn_0_920863N
8conv_1_kernel_regularizer_square_readvariableop_resource: 
identityЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpп
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp8conv_1_kernel_regularizer_square_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul
IdentityIdentity!Conv_1/kernel/Regularizer/mul:z:00^Conv_1/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp
н
_
C__inference_flatten_layer_call_and_return_conditional_losses_920720

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
и$
ѓ
B__inference_Conv_1_layer_call_and_return_conditional_losses_919141

inputsA
+conv1d_expanddims_1_readvariableop_resource: -
biasadd_readvariableop_resource: 
identityЂBiasAdd/ReadVariableOpЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂ"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
conv1d/ExpandDims/dim
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџШ2
conv1d/ExpandDimsИ
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dimЗ
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
conv1d/ExpandDims_1И
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:џџџџџџџџџА *
paddingVALID*
strides
2
conv1d
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:џџџџџџџџџА *
squeeze_dims

§џџџџџџџџ2
conv1d/Squeeze
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:џџџџџџџџџА 2	
BiasAdd]
TanhTanhBiasAdd:output:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
Tanhв
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulК
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:џџџџџџџџџА 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:џџџџџџџџџШ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs

Ц
/__inference_FCN_classifier_layer_call_fn_919769
convolutional_inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8 *S
fNRL
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_9197132
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
,
_output_shapes
:џџџџџџџџџШ
.
_user_specified_nameConvolutional_inputs

Ї
__inference_loss_fn_7_920940E
7dense_1_bias_regularizer_square_readvariableop_resource:@
identityЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpд
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp7dense_1_bias_regularizer_square_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul
IdentityIdentity Dense_1/bias/Regularizer/mul:z:0/^Dense_1/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp
Ѓ
L
0__inference_max_pooling1d_2_layer_call_fn_919106

inputs
identityп
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9191002
PartitionedCall
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ:e a
=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs


(__inference_Dense_2_layer_call_fn_920785

inputs
unknown:@
	unknown_0:
identityЂStatefulPartitionedCallѓ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9192992
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ@
 
_user_specified_nameinputs

a
C__inference_dropout_layer_call_and_return_conditional_losses_919153

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџX :S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs
Ћ

'__inference_Conv_3_layer_call_fn_920654

inputs
unknown:@ 
	unknown_0: 
identityЂStatefulPartitionedCallі
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9192252
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:џџџџџџџџџ  2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:џџџџџџџџџ&@: : 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
а$
ѓ
B__inference_Conv_2_layer_call_and_return_conditional_losses_920606

inputsA
+conv1d_expanddims_1_readvariableop_resource: @-
biasadd_readvariableop_resource:@
identityЂBiasAdd/ReadVariableOpЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂ"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
conv1d/ExpandDims/dim
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџX 2
conv1d/ExpandDimsИ
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dimЗ
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
conv1d/ExpandDims_1З
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@*
paddingVALID*
strides
2
conv1d
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@*
squeeze_dims

§џџџџџџџџ2
conv1d/Squeeze
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџL@2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
Tanhв
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulК
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:џџџџџџџџџL@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:џџџџџџџџџX : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs
Ъ
D
(__inference_dropout_layer_call_fn_920535

inputs
identityХ
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9191532
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџX :S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs

c
E__inference_dropout_1_layer_call_and_return_conditional_losses_919195

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ&@:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
П
й
C__inference_Dense_1_layer_call_and_return_conditional_losses_919270

inputs1
matmul_readvariableop_resource:	@-
biasadd_readvariableop_resource:@
identityЂBiasAdd/ReadVariableOpЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
TanhФ
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulМ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulё
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:џџџџџџџџџ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ь
d
E__inference_dropout_1_layer_call_and_return_conditional_losses_920633

inputs
identityc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/ShapeИ
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2
dropout/GreaterEqual/yТ
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџ&@2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ&@:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
я
Ѕ
__inference_loss_fn_5_920918D
6conv_3_bias_regularizer_square_readvariableop_resource: 
identityЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpб
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp6conv_3_bias_regularizer_square_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul
IdentityIdentityConv_3/bias/Regularizer/mul:z:0.^Conv_3/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp
Ѓ
L
0__inference_max_pooling1d_1_layer_call_fn_919091

inputs
identityп
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9190852
PartitionedCall
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ:e a
=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
Ъ
b
C__inference_dropout_layer_call_and_return_conditional_losses_919556

inputs
identityc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/ShapeИ
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџX *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2
dropout/GreaterEqual/yТ
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџX 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџX :S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs
Џ

'__inference_Conv_1_layer_call_fn_920502

inputs
unknown: 
	unknown_0: 
identityЂStatefulPartitionedCallї
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџА *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9191412
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:џџџџџџџџџА 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:џџџџџџџџџШ: : 22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
П
й
C__inference_Dense_1_layer_call_and_return_conditional_losses_920764

inputs1
matmul_readvariableop_resource:	@-
biasadd_readvariableop_resource:@
identityЂBiasAdd/ReadVariableOpЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
TanhФ
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulМ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulё
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:џџџџџџџџџ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ъ
b
C__inference_dropout_layer_call_and_return_conditional_losses_920557

inputs
identityc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/ShapeИ
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџX *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2
dropout/GreaterEqual/yТ
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџX 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџX :S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs

Ї
__inference_loss_fn_9_920962E
7dense_2_bias_regularizer_square_readvariableop_resource:
identityЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpд
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp7dense_2_bias_regularizer_square_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul
IdentityIdentity Dense_2/bias/Regularizer/mul:z:0/^Dense_2/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp

c
E__inference_dropout_1_layer_call_and_return_conditional_losses_920621

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ&@:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs

Ъ
__inference_loss_fn_10_920973X
Fdense_3_with_softmax_kernel_regularizer_square_readvariableop_resource:
identityЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpFdense_3_with_softmax_kernel_regularizer_square_readvariableop_resource*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulВ
IdentityIdentity/Dense_3_with_Softmax/kernel/Regularizer/mul:z:0>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp


(__inference_Dense_1_layer_call_fn_920741

inputs
unknown:	@
	unknown_0:@
identityЂStatefulPartitionedCallѓ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9192702
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:џџџџџџџџџ: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ћ

'__inference_Conv_2_layer_call_fn_920578

inputs
unknown: @
	unknown_0:@
identityЂStatefulPartitionedCallі
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџL@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9191832
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:џџџџџџџџџL@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:џџџџџџџџџX : : 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs

c
E__inference_dropout_2_layer_call_and_return_conditional_losses_919237

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
ќЌ
к

J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919882
convolutional_inputs#
conv_1_919772: 
conv_1_919774: #
conv_2_919779: @
conv_2_919781:@#
conv_3_919786:@ 
conv_3_919788: !
dense_1_919794:	@
dense_1_919796:@ 
dense_2_919799:@
dense_2_919801:-
dense_3_with_softmax_919804:)
dense_3_with_softmax_919806:
identityЂConv_1/StatefulPartitionedCallЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂConv_2/StatefulPartitionedCallЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂConv_3/StatefulPartitionedCallЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂDense_1/StatefulPartitionedCallЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂDense_2/StatefulPartitionedCallЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂ,Dense_3_with_Softmax/StatefulPartitionedCallЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp
Conv_1/StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsconv_1_919772conv_1_919774*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџА *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9191412 
Conv_1/StatefulPartitionedCall
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9190702
max_pooling1d/PartitionedCallѕ
dropout/PartitionedCallPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9191532
dropout/PartitionedCallЈ
Conv_2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0conv_2_919779conv_2_919781*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџL@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9191832 
Conv_2/StatefulPartitionedCall
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9190852!
max_pooling1d_1/PartitionedCall§
dropout_1/PartitionedCallPartitionedCall(max_pooling1d_1/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9191952
dropout_1/PartitionedCallЊ
Conv_3/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0conv_3_919786conv_3_919788*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9192252 
Conv_3/StatefulPartitionedCall
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9191002!
max_pooling1d_2/PartitionedCall§
dropout_2/PartitionedCallPartitionedCall(max_pooling1d_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9192372
dropout_2/PartitionedCallю
flatten/PartitionedCallPartitionedCall"dropout_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:џџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9192452
flatten/PartitionedCallЉ
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_919794dense_1_919796*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9192702!
Dense_1/StatefulPartitionedCallБ
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_919799dense_2_919801*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9192992!
Dense_2/StatefulPartitionedCallђ
,Dense_3_with_Softmax/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_with_softmax_919804dense_3_with_softmax_919806*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *Y
fTRR
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_9193282.
,Dense_3_with_Softmax/StatefulPartitionedCallД
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919772*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulЈ
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919774*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mulД
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919779*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulЈ
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919781*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mulД
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919786*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulЈ
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919788*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mulД
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919794*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulЋ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919796*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulГ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919799*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulЋ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919801*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulк
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919804*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulв
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919806*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЫ
IdentityIdentity5Dense_3_with_Softmax/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp-^Dense_3_with_Softmax/StatefulPartitionedCall<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2@
Conv_1/StatefulPartitionedCallConv_1/StatefulPartitionedCall2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2@
Conv_2/StatefulPartitionedCallConv_2/StatefulPartitionedCall2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2@
Conv_3/StatefulPartitionedCallConv_3/StatefulPartitionedCall2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2B
Dense_1/StatefulPartitionedCallDense_1/StatefulPartitionedCall2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2B
Dense_2/StatefulPartitionedCallDense_2/StatefulPartitionedCall2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2\
,Dense_3_with_Softmax/StatefulPartitionedCall,Dense_3_with_Softmax/StatefulPartitionedCall2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:b ^
,
_output_shapes
:џџџџџџџџџШ
.
_user_specified_nameConvolutional_inputs
Й
и
C__inference_Dense_2_layer_call_and_return_conditional_losses_919299

inputs0
matmul_readvariableop_resource:@-
biasadd_readvariableop_resource:
identityЂBiasAdd/ReadVariableOpЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
TanhУ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulМ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulё
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ@
 
_user_specified_nameinputs
Т
Б
__inference_loss_fn_4_920907N
8conv_3_kernel_regularizer_square_readvariableop_resource:@ 
identityЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpп
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp8conv_3_kernel_regularizer_square_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul
IdentityIdentity!Conv_3/kernel/Regularizer/mul:z:00^Conv_3/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp
ЪБ
Ф
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919995
convolutional_inputs#
conv_1_919885: 
conv_1_919887: #
conv_2_919892: @
conv_2_919894:@#
conv_3_919899:@ 
conv_3_919901: !
dense_1_919907:	@
dense_1_919909:@ 
dense_2_919912:@
dense_2_919914:-
dense_3_with_softmax_919917:)
dense_3_with_softmax_919919:
identityЂConv_1/StatefulPartitionedCallЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂConv_2/StatefulPartitionedCallЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂConv_3/StatefulPartitionedCallЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂDense_1/StatefulPartitionedCallЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂDense_2/StatefulPartitionedCallЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂ,Dense_3_with_Softmax/StatefulPartitionedCallЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpЂdropout/StatefulPartitionedCallЂ!dropout_1/StatefulPartitionedCallЂ!dropout_2/StatefulPartitionedCall
Conv_1/StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsconv_1_919885conv_1_919887*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџА *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9191412 
Conv_1/StatefulPartitionedCall
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9190702
max_pooling1d/PartitionedCall
dropout/StatefulPartitionedCallStatefulPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9195562!
dropout/StatefulPartitionedCallА
Conv_2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0conv_2_919892conv_2_919894*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџL@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9191832 
Conv_2/StatefulPartitionedCall
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9190852!
max_pooling1d_1/PartitionedCallЗ
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_1/PartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9195232#
!dropout_1/StatefulPartitionedCallВ
Conv_3/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0conv_3_919899conv_3_919901*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9192252 
Conv_3/StatefulPartitionedCall
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9191002!
max_pooling1d_2/PartitionedCallЙ
!dropout_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_2/PartitionedCall:output:0"^dropout_1/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9194902#
!dropout_2/StatefulPartitionedCallі
flatten/PartitionedCallPartitionedCall*dropout_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:џџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9192452
flatten/PartitionedCallЉ
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_919907dense_1_919909*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9192702!
Dense_1/StatefulPartitionedCallБ
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_919912dense_2_919914*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9192992!
Dense_2/StatefulPartitionedCallђ
,Dense_3_with_Softmax/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_with_softmax_919917dense_3_with_softmax_919919*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *Y
fTRR
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_9193282.
,Dense_3_with_Softmax/StatefulPartitionedCallД
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919885*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulЈ
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919887*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mulД
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919892*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulЈ
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919894*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mulД
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919899*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulЈ
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919901*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mulД
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919907*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulЋ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919909*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulГ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919912*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulЋ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919914*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulк
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919917*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulв
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919919*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЕ
IdentityIdentity5Dense_3_with_Softmax/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp-^Dense_3_with_Softmax/StatefulPartitionedCall<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp ^dropout/StatefulPartitionedCall"^dropout_1/StatefulPartitionedCall"^dropout_2/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2@
Conv_1/StatefulPartitionedCallConv_1/StatefulPartitionedCall2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2@
Conv_2/StatefulPartitionedCallConv_2/StatefulPartitionedCall2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2@
Conv_3/StatefulPartitionedCallConv_3/StatefulPartitionedCall2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2B
Dense_1/StatefulPartitionedCallDense_1/StatefulPartitionedCall2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2B
Dense_2/StatefulPartitionedCallDense_2/StatefulPartitionedCall2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2\
,Dense_3_with_Softmax/StatefulPartitionedCall,Dense_3_with_Softmax/StatefulPartitionedCall2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall2F
!dropout_2/StatefulPartitionedCall!dropout_2/StatefulPartitionedCall:b ^
,
_output_shapes
:џџџџџџџџџШ
.
_user_specified_nameConvolutional_inputs
 Б
Ж
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919713

inputs#
conv_1_919603: 
conv_1_919605: #
conv_2_919610: @
conv_2_919612:@#
conv_3_919617:@ 
conv_3_919619: !
dense_1_919625:	@
dense_1_919627:@ 
dense_2_919630:@
dense_2_919632:-
dense_3_with_softmax_919635:)
dense_3_with_softmax_919637:
identityЂConv_1/StatefulPartitionedCallЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂConv_2/StatefulPartitionedCallЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂConv_3/StatefulPartitionedCallЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂDense_1/StatefulPartitionedCallЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂDense_2/StatefulPartitionedCallЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂ,Dense_3_with_Softmax/StatefulPartitionedCallЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpЂdropout/StatefulPartitionedCallЂ!dropout_1/StatefulPartitionedCallЂ!dropout_2/StatefulPartitionedCall
Conv_1/StatefulPartitionedCallStatefulPartitionedCallinputsconv_1_919603conv_1_919605*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџА *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9191412 
Conv_1/StatefulPartitionedCall
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9190702
max_pooling1d/PartitionedCall
dropout/StatefulPartitionedCallStatefulPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9195562!
dropout/StatefulPartitionedCallА
Conv_2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0conv_2_919610conv_2_919612*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџL@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9191832 
Conv_2/StatefulPartitionedCall
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9190852!
max_pooling1d_1/PartitionedCallЗ
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_1/PartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9195232#
!dropout_1/StatefulPartitionedCallВ
Conv_3/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0conv_3_919617conv_3_919619*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9192252 
Conv_3/StatefulPartitionedCall
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9191002!
max_pooling1d_2/PartitionedCallЙ
!dropout_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_2/PartitionedCall:output:0"^dropout_1/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9194902#
!dropout_2/StatefulPartitionedCallі
flatten/PartitionedCallPartitionedCall*dropout_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:џџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9192452
flatten/PartitionedCallЉ
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_919625dense_1_919627*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9192702!
Dense_1/StatefulPartitionedCallБ
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_919630dense_2_919632*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9192992!
Dense_2/StatefulPartitionedCallђ
,Dense_3_with_Softmax/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_with_softmax_919635dense_3_with_softmax_919637*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *Y
fTRR
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_9193282.
,Dense_3_with_Softmax/StatefulPartitionedCallД
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919603*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulЈ
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919605*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mulД
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919610*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulЈ
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919612*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mulД
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919617*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulЈ
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919619*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mulД
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919625*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulЋ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919627*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulГ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919630*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulЋ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919632*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulк
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919635*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulв
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919637*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЕ
IdentityIdentity5Dense_3_with_Softmax/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp-^Dense_3_with_Softmax/StatefulPartitionedCall<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp ^dropout/StatefulPartitionedCall"^dropout_1/StatefulPartitionedCall"^dropout_2/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2@
Conv_1/StatefulPartitionedCallConv_1/StatefulPartitionedCall2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2@
Conv_2/StatefulPartitionedCallConv_2/StatefulPartitionedCall2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2@
Conv_3/StatefulPartitionedCallConv_3/StatefulPartitionedCall2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2B
Dense_1/StatefulPartitionedCallDense_1/StatefulPartitionedCall2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2B
Dense_2/StatefulPartitionedCallDense_2/StatefulPartitionedCall2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2\
,Dense_3_with_Softmax/StatefulPartitionedCall,Dense_3_with_Softmax/StatefulPartitionedCall2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall2F
!dropout_2/StatefulPartitionedCall!dropout_2/StatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
Э
А
__inference_loss_fn_6_920929L
9dense_1_kernel_regularizer_square_readvariableop_resource:	@
identityЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpп
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_1_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul
IdentityIdentity"Dense_1/kernel/Regularizer/mul:z:01^Dense_1/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp
 §

J__inference_FCN_classifier_layer_call_and_return_conditional_losses_920481

inputsH
2conv_1_conv1d_expanddims_1_readvariableop_resource: 4
&conv_1_biasadd_readvariableop_resource: H
2conv_2_conv1d_expanddims_1_readvariableop_resource: @4
&conv_2_biasadd_readvariableop_resource:@H
2conv_3_conv1d_expanddims_1_readvariableop_resource:@ 4
&conv_3_biasadd_readvariableop_resource: 9
&dense_1_matmul_readvariableop_resource:	@5
'dense_1_biasadd_readvariableop_resource:@8
&dense_2_matmul_readvariableop_resource:@5
'dense_2_biasadd_readvariableop_resource:E
3dense_3_with_softmax_matmul_readvariableop_resource:B
4dense_3_with_softmax_biasadd_readvariableop_resource:
identityЂConv_1/BiasAdd/ReadVariableOpЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ)Conv_1/conv1d/ExpandDims_1/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂConv_2/BiasAdd/ReadVariableOpЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ)Conv_2/conv1d/ExpandDims_1/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂConv_3/BiasAdd/ReadVariableOpЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ)Conv_3/conv1d/ExpandDims_1/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂDense_1/BiasAdd/ReadVariableOpЂDense_1/MatMul/ReadVariableOpЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂDense_2/BiasAdd/ReadVariableOpЂDense_2/MatMul/ReadVariableOpЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂ+Dense_3_with_Softmax/BiasAdd/ReadVariableOpЂ*Dense_3_with_Softmax/MatMul/ReadVariableOpЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp
Conv_1/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
Conv_1/conv1d/ExpandDims/dimЌ
Conv_1/conv1d/ExpandDims
ExpandDimsinputs%Conv_1/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџШ2
Conv_1/conv1d/ExpandDimsЭ
)Conv_1/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02+
)Conv_1/conv1d/ExpandDims_1/ReadVariableOp
Conv_1/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_1/conv1d/ExpandDims_1/dimг
Conv_1/conv1d/ExpandDims_1
ExpandDims1Conv_1/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_1/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
Conv_1/conv1d/ExpandDims_1д
Conv_1/conv1dConv2D!Conv_1/conv1d/ExpandDims:output:0#Conv_1/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:џџџџџџџџџА *
paddingVALID*
strides
2
Conv_1/conv1dЈ
Conv_1/conv1d/SqueezeSqueezeConv_1/conv1d:output:0*
T0*,
_output_shapes
:џџџџџџџџџА *
squeeze_dims

§џџџџџџџџ2
Conv_1/conv1d/SqueezeЁ
Conv_1/BiasAdd/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_1/BiasAdd/ReadVariableOpЉ
Conv_1/BiasAddBiasAddConv_1/conv1d/Squeeze:output:0%Conv_1/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
Conv_1/BiasAddr
Conv_1/TanhTanhConv_1/BiasAdd:output:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
Conv_1/Tanh~
max_pooling1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
max_pooling1d/ExpandDims/dimЕ
max_pooling1d/ExpandDims
ExpandDimsConv_1/Tanh:y:0%max_pooling1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџА 2
max_pooling1d/ExpandDimsЩ
max_pooling1d/MaxPoolMaxPool!max_pooling1d/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџX *
ksize
*
paddingVALID*
strides
2
max_pooling1d/MaxPoolІ
max_pooling1d/SqueezeSqueezemax_pooling1d/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџX *
squeeze_dims
2
max_pooling1d/Squeezes
dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/dropout/ConstЇ
dropout/dropout/MulMulmax_pooling1d/Squeeze:output:0dropout/dropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/dropout/Mul|
dropout/dropout/ShapeShapemax_pooling1d/Squeeze:output:0*
T0*
_output_shapes
:2
dropout/dropout/Shapeа
,dropout/dropout/random_uniform/RandomUniformRandomUniformdropout/dropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџX *
dtype02.
,dropout/dropout/random_uniform/RandomUniform
dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2 
dropout/dropout/GreaterEqual/yт
dropout/dropout/GreaterEqualGreaterEqual5dropout/dropout/random_uniform/RandomUniform:output:0'dropout/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/dropout/GreaterEqual
dropout/dropout/CastCast dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџX 2
dropout/dropout/Cast
dropout/dropout/Mul_1Muldropout/dropout/Mul:z:0dropout/dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/dropout/Mul_1
Conv_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
Conv_2/conv1d/ExpandDims/dimО
Conv_2/conv1d/ExpandDims
ExpandDimsdropout/dropout/Mul_1:z:0%Conv_2/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџX 2
Conv_2/conv1d/ExpandDimsЭ
)Conv_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02+
)Conv_2/conv1d/ExpandDims_1/ReadVariableOp
Conv_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_2/conv1d/ExpandDims_1/dimг
Conv_2/conv1d/ExpandDims_1
ExpandDims1Conv_2/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_2/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
Conv_2/conv1d/ExpandDims_1г
Conv_2/conv1dConv2D!Conv_2/conv1d/ExpandDims:output:0#Conv_2/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@*
paddingVALID*
strides
2
Conv_2/conv1dЇ
Conv_2/conv1d/SqueezeSqueezeConv_2/conv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@*
squeeze_dims

§џџџџџџџџ2
Conv_2/conv1d/SqueezeЁ
Conv_2/BiasAdd/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
Conv_2/BiasAdd/ReadVariableOpЈ
Conv_2/BiasAddBiasAddConv_2/conv1d/Squeeze:output:0%Conv_2/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
Conv_2/BiasAddq
Conv_2/TanhTanhConv_2/BiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
Conv_2/Tanh
max_pooling1d_1/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_1/ExpandDims/dimК
max_pooling1d_1/ExpandDims
ExpandDimsConv_2/Tanh:y:0'max_pooling1d_1/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@2
max_pooling1d_1/ExpandDimsЯ
max_pooling1d_1/MaxPoolMaxPool#max_pooling1d_1/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџ&@*
ksize
*
paddingVALID*
strides
2
max_pooling1d_1/MaxPoolЌ
max_pooling1d_1/SqueezeSqueeze max_pooling1d_1/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@*
squeeze_dims
2
max_pooling1d_1/Squeezew
dropout_1/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout_1/dropout/ConstЏ
dropout_1/dropout/MulMul max_pooling1d_1/Squeeze:output:0 dropout_1/dropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout_1/dropout/Mul
dropout_1/dropout/ShapeShape max_pooling1d_1/Squeeze:output:0*
T0*
_output_shapes
:2
dropout_1/dropout/Shapeж
.dropout_1/dropout/random_uniform/RandomUniformRandomUniform dropout_1/dropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@*
dtype020
.dropout_1/dropout/random_uniform/RandomUniform
 dropout_1/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2"
 dropout_1/dropout/GreaterEqual/yъ
dropout_1/dropout/GreaterEqualGreaterEqual7dropout_1/dropout/random_uniform/RandomUniform:output:0)dropout_1/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2 
dropout_1/dropout/GreaterEqualЁ
dropout_1/dropout/CastCast"dropout_1/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџ&@2
dropout_1/dropout/CastІ
dropout_1/dropout/Mul_1Muldropout_1/dropout/Mul:z:0dropout_1/dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout_1/dropout/Mul_1
Conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
Conv_3/conv1d/ExpandDims/dimР
Conv_3/conv1d/ExpandDims
ExpandDimsdropout_1/dropout/Mul_1:z:0%Conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ&@2
Conv_3/conv1d/ExpandDimsЭ
)Conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02+
)Conv_3/conv1d/ExpandDims_1/ReadVariableOp
Conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_3/conv1d/ExpandDims_1/dimг
Conv_3/conv1d/ExpandDims_1
ExpandDims1Conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
Conv_3/conv1d/ExpandDims_1г
Conv_3/conv1dConv2D!Conv_3/conv1d/ExpandDims:output:0#Conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  *
paddingVALID*
strides
2
Conv_3/conv1dЇ
Conv_3/conv1d/SqueezeSqueezeConv_3/conv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  *
squeeze_dims

§џџџџџџџџ2
Conv_3/conv1d/SqueezeЁ
Conv_3/BiasAdd/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_3/BiasAdd/ReadVariableOpЈ
Conv_3/BiasAddBiasAddConv_3/conv1d/Squeeze:output:0%Conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
Conv_3/BiasAddq
Conv_3/TanhTanhConv_3/BiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
Conv_3/Tanh
max_pooling1d_2/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_2/ExpandDims/dimК
max_pooling1d_2/ExpandDims
ExpandDimsConv_3/Tanh:y:0'max_pooling1d_2/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  2
max_pooling1d_2/ExpandDimsЯ
max_pooling1d_2/MaxPoolMaxPool#max_pooling1d_2/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџ *
ksize
*
paddingVALID*
strides
2
max_pooling1d_2/MaxPoolЌ
max_pooling1d_2/SqueezeSqueeze max_pooling1d_2/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџ *
squeeze_dims
2
max_pooling1d_2/Squeezew
dropout_2/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout_2/dropout/ConstЏ
dropout_2/dropout/MulMul max_pooling1d_2/Squeeze:output:0 dropout_2/dropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout_2/dropout/Mul
dropout_2/dropout/ShapeShape max_pooling1d_2/Squeeze:output:0*
T0*
_output_shapes
:2
dropout_2/dropout/Shapeж
.dropout_2/dropout/random_uniform/RandomUniformRandomUniform dropout_2/dropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџ *
dtype020
.dropout_2/dropout/random_uniform/RandomUniform
 dropout_2/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2"
 dropout_2/dropout/GreaterEqual/yъ
dropout_2/dropout/GreaterEqualGreaterEqual7dropout_2/dropout/random_uniform/RandomUniform:output:0)dropout_2/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2 
dropout_2/dropout/GreaterEqualЁ
dropout_2/dropout/CastCast"dropout_2/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџ 2
dropout_2/dropout/CastІ
dropout_2/dropout/Mul_1Muldropout_2/dropout/Mul:z:0dropout_2/dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout_2/dropout/Mul_1o
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
flatten/Const
flatten/ReshapeReshapedropout_2/dropout/Mul_1:z:0flatten/Const:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2
flatten/ReshapeІ
Dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	@*
dtype02
Dense_1/MatMul/ReadVariableOp
Dense_1/MatMulMatMulflatten/Reshape:output:0%Dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
Dense_1/MatMulЄ
Dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02 
Dense_1/BiasAdd/ReadVariableOpЁ
Dense_1/BiasAddBiasAddDense_1/MatMul:product:0&Dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
Dense_1/BiasAddp
Dense_1/TanhTanhDense_1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
Dense_1/TanhЅ
Dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype02
Dense_2/MatMul/ReadVariableOp
Dense_2/MatMulMatMulDense_1/Tanh:y:0%Dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_2/MatMulЄ
Dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
Dense_2/BiasAdd/ReadVariableOpЁ
Dense_2/BiasAddBiasAddDense_2/MatMul:product:0&Dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_2/BiasAddp
Dense_2/TanhTanhDense_2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_2/TanhЬ
*Dense_3_with_Softmax/MatMul/ReadVariableOpReadVariableOp3dense_3_with_softmax_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*Dense_3_with_Softmax/MatMul/ReadVariableOpМ
Dense_3_with_Softmax/MatMulMatMulDense_2/Tanh:y:02Dense_3_with_Softmax/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_3_with_Softmax/MatMulЫ
+Dense_3_with_Softmax/BiasAdd/ReadVariableOpReadVariableOp4dense_3_with_softmax_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+Dense_3_with_Softmax/BiasAdd/ReadVariableOpе
Dense_3_with_Softmax/BiasAddBiasAdd%Dense_3_with_Softmax/MatMul:product:03Dense_3_with_Softmax/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_3_with_Softmax/BiasAdd 
Dense_3_with_Softmax/SoftmaxSoftmax%Dense_3_with_Softmax/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_3_with_Softmax/Softmaxй
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulС
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mulй
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulС
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mulй
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulС
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mulЬ
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulФ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulЫ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulФ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulђ
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOp3dense_3_with_softmax_matmul_readvariableop_resource*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulы
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOp4dense_3_with_softmax_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЇ	
IdentityIdentity&Dense_3_with_Softmax/Softmax:softmax:0^Conv_1/BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp*^Conv_1/conv1d/ExpandDims_1/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp*^Conv_2/conv1d/ExpandDims_1/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp*^Conv_3/conv1d/ExpandDims_1/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp^Dense_1/BiasAdd/ReadVariableOp^Dense_1/MatMul/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^Dense_2/BiasAdd/ReadVariableOp^Dense_2/MatMul/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp,^Dense_3_with_Softmax/BiasAdd/ReadVariableOp+^Dense_3_with_Softmax/MatMul/ReadVariableOp<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2>
Conv_1/BiasAdd/ReadVariableOpConv_1/BiasAdd/ReadVariableOp2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2V
)Conv_1/conv1d/ExpandDims_1/ReadVariableOp)Conv_1/conv1d/ExpandDims_1/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2>
Conv_2/BiasAdd/ReadVariableOpConv_2/BiasAdd/ReadVariableOp2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2V
)Conv_2/conv1d/ExpandDims_1/ReadVariableOp)Conv_2/conv1d/ExpandDims_1/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2>
Conv_3/BiasAdd/ReadVariableOpConv_3/BiasAdd/ReadVariableOp2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2V
)Conv_3/conv1d/ExpandDims_1/ReadVariableOp)Conv_3/conv1d/ExpandDims_1/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2@
Dense_1/BiasAdd/ReadVariableOpDense_1/BiasAdd/ReadVariableOp2>
Dense_1/MatMul/ReadVariableOpDense_1/MatMul/ReadVariableOp2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2@
Dense_2/BiasAdd/ReadVariableOpDense_2/BiasAdd/ReadVariableOp2>
Dense_2/MatMul/ReadVariableOpDense_2/MatMul/ReadVariableOp2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2Z
+Dense_3_with_Softmax/BiasAdd/ReadVariableOp+Dense_3_with_Softmax/BiasAdd/ReadVariableOp2X
*Dense_3_with_Softmax/MatMul/ReadVariableOp*Dense_3_with_Softmax/MatMul/ReadVariableOp2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
я
Ѕ
__inference_loss_fn_1_920874D
6conv_1_bias_regularizer_square_readvariableop_resource: 
identityЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpб
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp6conv_1_bias_regularizer_square_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul
IdentityIdentityConv_1/bias/Regularizer/mul:z:0.^Conv_1/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp
О
х
!__inference__wrapped_model_919061
convolutional_inputsW
Afcn_classifier_conv_1_conv1d_expanddims_1_readvariableop_resource: C
5fcn_classifier_conv_1_biasadd_readvariableop_resource: W
Afcn_classifier_conv_2_conv1d_expanddims_1_readvariableop_resource: @C
5fcn_classifier_conv_2_biasadd_readvariableop_resource:@W
Afcn_classifier_conv_3_conv1d_expanddims_1_readvariableop_resource:@ C
5fcn_classifier_conv_3_biasadd_readvariableop_resource: H
5fcn_classifier_dense_1_matmul_readvariableop_resource:	@D
6fcn_classifier_dense_1_biasadd_readvariableop_resource:@G
5fcn_classifier_dense_2_matmul_readvariableop_resource:@D
6fcn_classifier_dense_2_biasadd_readvariableop_resource:T
Bfcn_classifier_dense_3_with_softmax_matmul_readvariableop_resource:Q
Cfcn_classifier_dense_3_with_softmax_biasadd_readvariableop_resource:
identityЂ,FCN_classifier/Conv_1/BiasAdd/ReadVariableOpЂ8FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOpЂ,FCN_classifier/Conv_2/BiasAdd/ReadVariableOpЂ8FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOpЂ,FCN_classifier/Conv_3/BiasAdd/ReadVariableOpЂ8FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOpЂ-FCN_classifier/Dense_1/BiasAdd/ReadVariableOpЂ,FCN_classifier/Dense_1/MatMul/ReadVariableOpЂ-FCN_classifier/Dense_2/BiasAdd/ReadVariableOpЂ,FCN_classifier/Dense_2/MatMul/ReadVariableOpЂ:FCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOpЂ9FCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOpЅ
+FCN_classifier/Conv_1/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2-
+FCN_classifier/Conv_1/conv1d/ExpandDims/dimч
'FCN_classifier/Conv_1/conv1d/ExpandDims
ExpandDimsconvolutional_inputs4FCN_classifier/Conv_1/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџШ2)
'FCN_classifier/Conv_1/conv1d/ExpandDimsњ
8FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOpReadVariableOpAfcn_classifier_conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02:
8FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOp 
-FCN_classifier/Conv_1/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2/
-FCN_classifier/Conv_1/conv1d/ExpandDims_1/dim
)FCN_classifier/Conv_1/conv1d/ExpandDims_1
ExpandDims@FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOp:value:06FCN_classifier/Conv_1/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2+
)FCN_classifier/Conv_1/conv1d/ExpandDims_1
FCN_classifier/Conv_1/conv1dConv2D0FCN_classifier/Conv_1/conv1d/ExpandDims:output:02FCN_classifier/Conv_1/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:џџџџџџџџџА *
paddingVALID*
strides
2
FCN_classifier/Conv_1/conv1dе
$FCN_classifier/Conv_1/conv1d/SqueezeSqueeze%FCN_classifier/Conv_1/conv1d:output:0*
T0*,
_output_shapes
:џџџџџџџџџА *
squeeze_dims

§џџџџџџџџ2&
$FCN_classifier/Conv_1/conv1d/SqueezeЮ
,FCN_classifier/Conv_1/BiasAdd/ReadVariableOpReadVariableOp5fcn_classifier_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02.
,FCN_classifier/Conv_1/BiasAdd/ReadVariableOpх
FCN_classifier/Conv_1/BiasAddBiasAdd-FCN_classifier/Conv_1/conv1d/Squeeze:output:04FCN_classifier/Conv_1/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
FCN_classifier/Conv_1/BiasAdd
FCN_classifier/Conv_1/TanhTanh&FCN_classifier/Conv_1/BiasAdd:output:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
FCN_classifier/Conv_1/Tanh
+FCN_classifier/max_pooling1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2-
+FCN_classifier/max_pooling1d/ExpandDims/dimё
'FCN_classifier/max_pooling1d/ExpandDims
ExpandDimsFCN_classifier/Conv_1/Tanh:y:04FCN_classifier/max_pooling1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџА 2)
'FCN_classifier/max_pooling1d/ExpandDimsі
$FCN_classifier/max_pooling1d/MaxPoolMaxPool0FCN_classifier/max_pooling1d/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџX *
ksize
*
paddingVALID*
strides
2&
$FCN_classifier/max_pooling1d/MaxPoolг
$FCN_classifier/max_pooling1d/SqueezeSqueeze-FCN_classifier/max_pooling1d/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџX *
squeeze_dims
2&
$FCN_classifier/max_pooling1d/SqueezeГ
FCN_classifier/dropout/IdentityIdentity-FCN_classifier/max_pooling1d/Squeeze:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2!
FCN_classifier/dropout/IdentityЅ
+FCN_classifier/Conv_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2-
+FCN_classifier/Conv_2/conv1d/ExpandDims/dimњ
'FCN_classifier/Conv_2/conv1d/ExpandDims
ExpandDims(FCN_classifier/dropout/Identity:output:04FCN_classifier/Conv_2/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџX 2)
'FCN_classifier/Conv_2/conv1d/ExpandDimsњ
8FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOpAfcn_classifier_conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02:
8FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOp 
-FCN_classifier/Conv_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2/
-FCN_classifier/Conv_2/conv1d/ExpandDims_1/dim
)FCN_classifier/Conv_2/conv1d/ExpandDims_1
ExpandDims@FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOp:value:06FCN_classifier/Conv_2/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2+
)FCN_classifier/Conv_2/conv1d/ExpandDims_1
FCN_classifier/Conv_2/conv1dConv2D0FCN_classifier/Conv_2/conv1d/ExpandDims:output:02FCN_classifier/Conv_2/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@*
paddingVALID*
strides
2
FCN_classifier/Conv_2/conv1dд
$FCN_classifier/Conv_2/conv1d/SqueezeSqueeze%FCN_classifier/Conv_2/conv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@*
squeeze_dims

§џџџџџџџџ2&
$FCN_classifier/Conv_2/conv1d/SqueezeЮ
,FCN_classifier/Conv_2/BiasAdd/ReadVariableOpReadVariableOp5fcn_classifier_conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02.
,FCN_classifier/Conv_2/BiasAdd/ReadVariableOpф
FCN_classifier/Conv_2/BiasAddBiasAdd-FCN_classifier/Conv_2/conv1d/Squeeze:output:04FCN_classifier/Conv_2/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
FCN_classifier/Conv_2/BiasAdd
FCN_classifier/Conv_2/TanhTanh&FCN_classifier/Conv_2/BiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
FCN_classifier/Conv_2/Tanh 
-FCN_classifier/max_pooling1d_1/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2/
-FCN_classifier/max_pooling1d_1/ExpandDims/dimі
)FCN_classifier/max_pooling1d_1/ExpandDims
ExpandDimsFCN_classifier/Conv_2/Tanh:y:06FCN_classifier/max_pooling1d_1/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@2+
)FCN_classifier/max_pooling1d_1/ExpandDimsќ
&FCN_classifier/max_pooling1d_1/MaxPoolMaxPool2FCN_classifier/max_pooling1d_1/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџ&@*
ksize
*
paddingVALID*
strides
2(
&FCN_classifier/max_pooling1d_1/MaxPoolй
&FCN_classifier/max_pooling1d_1/SqueezeSqueeze/FCN_classifier/max_pooling1d_1/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@*
squeeze_dims
2(
&FCN_classifier/max_pooling1d_1/SqueezeЙ
!FCN_classifier/dropout_1/IdentityIdentity/FCN_classifier/max_pooling1d_1/Squeeze:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2#
!FCN_classifier/dropout_1/IdentityЅ
+FCN_classifier/Conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2-
+FCN_classifier/Conv_3/conv1d/ExpandDims/dimќ
'FCN_classifier/Conv_3/conv1d/ExpandDims
ExpandDims*FCN_classifier/dropout_1/Identity:output:04FCN_classifier/Conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ&@2)
'FCN_classifier/Conv_3/conv1d/ExpandDimsњ
8FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOpAfcn_classifier_conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02:
8FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOp 
-FCN_classifier/Conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2/
-FCN_classifier/Conv_3/conv1d/ExpandDims_1/dim
)FCN_classifier/Conv_3/conv1d/ExpandDims_1
ExpandDims@FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:06FCN_classifier/Conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2+
)FCN_classifier/Conv_3/conv1d/ExpandDims_1
FCN_classifier/Conv_3/conv1dConv2D0FCN_classifier/Conv_3/conv1d/ExpandDims:output:02FCN_classifier/Conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  *
paddingVALID*
strides
2
FCN_classifier/Conv_3/conv1dд
$FCN_classifier/Conv_3/conv1d/SqueezeSqueeze%FCN_classifier/Conv_3/conv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  *
squeeze_dims

§џџџџџџџџ2&
$FCN_classifier/Conv_3/conv1d/SqueezeЮ
,FCN_classifier/Conv_3/BiasAdd/ReadVariableOpReadVariableOp5fcn_classifier_conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02.
,FCN_classifier/Conv_3/BiasAdd/ReadVariableOpф
FCN_classifier/Conv_3/BiasAddBiasAdd-FCN_classifier/Conv_3/conv1d/Squeeze:output:04FCN_classifier/Conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
FCN_classifier/Conv_3/BiasAdd
FCN_classifier/Conv_3/TanhTanh&FCN_classifier/Conv_3/BiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
FCN_classifier/Conv_3/Tanh 
-FCN_classifier/max_pooling1d_2/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2/
-FCN_classifier/max_pooling1d_2/ExpandDims/dimі
)FCN_classifier/max_pooling1d_2/ExpandDims
ExpandDimsFCN_classifier/Conv_3/Tanh:y:06FCN_classifier/max_pooling1d_2/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  2+
)FCN_classifier/max_pooling1d_2/ExpandDimsќ
&FCN_classifier/max_pooling1d_2/MaxPoolMaxPool2FCN_classifier/max_pooling1d_2/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџ *
ksize
*
paddingVALID*
strides
2(
&FCN_classifier/max_pooling1d_2/MaxPoolй
&FCN_classifier/max_pooling1d_2/SqueezeSqueeze/FCN_classifier/max_pooling1d_2/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџ *
squeeze_dims
2(
&FCN_classifier/max_pooling1d_2/SqueezeЙ
!FCN_classifier/dropout_2/IdentityIdentity/FCN_classifier/max_pooling1d_2/Squeeze:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2#
!FCN_classifier/dropout_2/Identity
FCN_classifier/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
FCN_classifier/flatten/Constб
FCN_classifier/flatten/ReshapeReshape*FCN_classifier/dropout_2/Identity:output:0%FCN_classifier/flatten/Const:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2 
FCN_classifier/flatten/Reshapeг
,FCN_classifier/Dense_1/MatMul/ReadVariableOpReadVariableOp5fcn_classifier_dense_1_matmul_readvariableop_resource*
_output_shapes
:	@*
dtype02.
,FCN_classifier/Dense_1/MatMul/ReadVariableOpй
FCN_classifier/Dense_1/MatMulMatMul'FCN_classifier/flatten/Reshape:output:04FCN_classifier/Dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
FCN_classifier/Dense_1/MatMulб
-FCN_classifier/Dense_1/BiasAdd/ReadVariableOpReadVariableOp6fcn_classifier_dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-FCN_classifier/Dense_1/BiasAdd/ReadVariableOpн
FCN_classifier/Dense_1/BiasAddBiasAdd'FCN_classifier/Dense_1/MatMul:product:05FCN_classifier/Dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2 
FCN_classifier/Dense_1/BiasAdd
FCN_classifier/Dense_1/TanhTanh'FCN_classifier/Dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
FCN_classifier/Dense_1/Tanhв
,FCN_classifier/Dense_2/MatMul/ReadVariableOpReadVariableOp5fcn_classifier_dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype02.
,FCN_classifier/Dense_2/MatMul/ReadVariableOpб
FCN_classifier/Dense_2/MatMulMatMulFCN_classifier/Dense_1/Tanh:y:04FCN_classifier/Dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
FCN_classifier/Dense_2/MatMulб
-FCN_classifier/Dense_2/BiasAdd/ReadVariableOpReadVariableOp6fcn_classifier_dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02/
-FCN_classifier/Dense_2/BiasAdd/ReadVariableOpн
FCN_classifier/Dense_2/BiasAddBiasAdd'FCN_classifier/Dense_2/MatMul:product:05FCN_classifier/Dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
FCN_classifier/Dense_2/BiasAdd
FCN_classifier/Dense_2/TanhTanh'FCN_classifier/Dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
FCN_classifier/Dense_2/Tanhљ
9FCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOpReadVariableOpBfcn_classifier_dense_3_with_softmax_matmul_readvariableop_resource*
_output_shapes

:*
dtype02;
9FCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOpј
*FCN_classifier/Dense_3_with_Softmax/MatMulMatMulFCN_classifier/Dense_2/Tanh:y:0AFCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2,
*FCN_classifier/Dense_3_with_Softmax/MatMulј
:FCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOpReadVariableOpCfcn_classifier_dense_3_with_softmax_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02<
:FCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOp
+FCN_classifier/Dense_3_with_Softmax/BiasAddBiasAdd4FCN_classifier/Dense_3_with_Softmax/MatMul:product:0BFCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+FCN_classifier/Dense_3_with_Softmax/BiasAddЭ
+FCN_classifier/Dense_3_with_Softmax/SoftmaxSoftmax4FCN_classifier/Dense_3_with_Softmax/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+FCN_classifier/Dense_3_with_Softmax/Softmaxў
IdentityIdentity5FCN_classifier/Dense_3_with_Softmax/Softmax:softmax:0-^FCN_classifier/Conv_1/BiasAdd/ReadVariableOp9^FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOp-^FCN_classifier/Conv_2/BiasAdd/ReadVariableOp9^FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOp-^FCN_classifier/Conv_3/BiasAdd/ReadVariableOp9^FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOp.^FCN_classifier/Dense_1/BiasAdd/ReadVariableOp-^FCN_classifier/Dense_1/MatMul/ReadVariableOp.^FCN_classifier/Dense_2/BiasAdd/ReadVariableOp-^FCN_classifier/Dense_2/MatMul/ReadVariableOp;^FCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOp:^FCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2\
,FCN_classifier/Conv_1/BiasAdd/ReadVariableOp,FCN_classifier/Conv_1/BiasAdd/ReadVariableOp2t
8FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOp8FCN_classifier/Conv_1/conv1d/ExpandDims_1/ReadVariableOp2\
,FCN_classifier/Conv_2/BiasAdd/ReadVariableOp,FCN_classifier/Conv_2/BiasAdd/ReadVariableOp2t
8FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOp8FCN_classifier/Conv_2/conv1d/ExpandDims_1/ReadVariableOp2\
,FCN_classifier/Conv_3/BiasAdd/ReadVariableOp,FCN_classifier/Conv_3/BiasAdd/ReadVariableOp2t
8FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOp8FCN_classifier/Conv_3/conv1d/ExpandDims_1/ReadVariableOp2^
-FCN_classifier/Dense_1/BiasAdd/ReadVariableOp-FCN_classifier/Dense_1/BiasAdd/ReadVariableOp2\
,FCN_classifier/Dense_1/MatMul/ReadVariableOp,FCN_classifier/Dense_1/MatMul/ReadVariableOp2^
-FCN_classifier/Dense_2/BiasAdd/ReadVariableOp-FCN_classifier/Dense_2/BiasAdd/ReadVariableOp2\
,FCN_classifier/Dense_2/MatMul/ReadVariableOp,FCN_classifier/Dense_2/MatMul/ReadVariableOp2x
:FCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOp:FCN_classifier/Dense_3_with_Softmax/BiasAdd/ReadVariableOp2v
9FCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOp9FCN_classifier/Dense_3_with_Softmax/MatMul/ReadVariableOp:b ^
,
_output_shapes
:џџџџџџџџџШ
.
_user_specified_nameConvolutional_inputs
а$
ѓ
B__inference_Conv_3_layer_call_and_return_conditional_losses_919225

inputsA
+conv1d_expanddims_1_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identityЂBiasAdd/ReadVariableOpЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂ"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
conv1d/ExpandDims/dim
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ&@2
conv1d/ExpandDimsИ
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dimЗ
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
conv1d/ExpandDims_1З
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  *
paddingVALID*
strides
2
conv1d
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  *
squeeze_dims

§џџџџџџџџ2
conv1d/Squeeze
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџ  2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
Tanhв
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulК
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:џџџџџџџџџ  2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:џџџџџџџџџ&@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs

Ц
/__inference_FCN_classifier_layer_call_fn_919434
convolutional_inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8 *S
fNRL
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_9194072
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
,
_output_shapes
:џџџџџџџџџШ
.
_user_specified_nameConvolutional_inputs
Г
Ђ
5__inference_Dense_3_with_Softmax_layer_call_fn_920829

inputs
unknown:
	unknown_0:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *Y
fTRR
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_9193282
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
њ!
џ
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_919328

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identityЂBiasAdd/ReadVariableOpЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAdda
SoftmaxSoftmaxBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Softmaxн
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulж
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mul
IdentityIdentitySoftmax:softmax:0^BiasAdd/ReadVariableOp<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Р
Т
__inference_loss_fn_11_920984R
Ddense_3_with_softmax_bias_regularizer_square_readvariableop_resource:
identityЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpћ
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpDdense_3_with_softmax_bias_regularizer_square_readvariableop_resource*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЎ
IdentityIdentity-Dense_3_with_Softmax/bias/Regularizer/mul:z:0<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp
Й
и
C__inference_Dense_2_layer_call_and_return_conditional_losses_920808

inputs0
matmul_readvariableop_resource:@-
biasadd_readvariableop_resource:
identityЂBiasAdd/ReadVariableOpЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂMatMul/ReadVariableOp
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
TanhУ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulМ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulё
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ@
 
_user_specified_nameinputs
Ъ
Џ
__inference_loss_fn_8_920951K
9dense_2_kernel_regularizer_square_readvariableop_resource:@
identityЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpо
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_2_kernel_regularizer_square_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul
IdentityIdentity"Dense_2/kernel/Regularizer/mul:z:01^Dense_2/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp

a
C__inference_dropout_layer_call_and_return_conditional_losses_920545

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџX :S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs
Ю
F
*__inference_dropout_1_layer_call_fn_920611

inputs
identityЧ
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9191952
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ&@:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
а$
ѓ
B__inference_Conv_2_layer_call_and_return_conditional_losses_919183

inputsA
+conv1d_expanddims_1_readvariableop_resource: @-
biasadd_readvariableop_resource:@
identityЂBiasAdd/ReadVariableOpЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂ"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
conv1d/ExpandDims/dim
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџX 2
conv1d/ExpandDimsИ
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dimЗ
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
conv1d/ExpandDims_1З
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@*
paddingVALID*
strides
2
conv1d
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@*
squeeze_dims

§џџџџџџџџ2
conv1d/Squeeze
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџL@2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
Tanhв
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulК
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:џџџџџџџџџL@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:џџџџџџџџџX : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs
с

J__inference_FCN_classifier_layer_call_and_return_conditional_losses_920310

inputsH
2conv_1_conv1d_expanddims_1_readvariableop_resource: 4
&conv_1_biasadd_readvariableop_resource: H
2conv_2_conv1d_expanddims_1_readvariableop_resource: @4
&conv_2_biasadd_readvariableop_resource:@H
2conv_3_conv1d_expanddims_1_readvariableop_resource:@ 4
&conv_3_biasadd_readvariableop_resource: 9
&dense_1_matmul_readvariableop_resource:	@5
'dense_1_biasadd_readvariableop_resource:@8
&dense_2_matmul_readvariableop_resource:@5
'dense_2_biasadd_readvariableop_resource:E
3dense_3_with_softmax_matmul_readvariableop_resource:B
4dense_3_with_softmax_biasadd_readvariableop_resource:
identityЂConv_1/BiasAdd/ReadVariableOpЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ)Conv_1/conv1d/ExpandDims_1/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂConv_2/BiasAdd/ReadVariableOpЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ)Conv_2/conv1d/ExpandDims_1/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂConv_3/BiasAdd/ReadVariableOpЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ)Conv_3/conv1d/ExpandDims_1/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂDense_1/BiasAdd/ReadVariableOpЂDense_1/MatMul/ReadVariableOpЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂDense_2/BiasAdd/ReadVariableOpЂDense_2/MatMul/ReadVariableOpЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂ+Dense_3_with_Softmax/BiasAdd/ReadVariableOpЂ*Dense_3_with_Softmax/MatMul/ReadVariableOpЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp
Conv_1/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
Conv_1/conv1d/ExpandDims/dimЌ
Conv_1/conv1d/ExpandDims
ExpandDimsinputs%Conv_1/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџШ2
Conv_1/conv1d/ExpandDimsЭ
)Conv_1/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02+
)Conv_1/conv1d/ExpandDims_1/ReadVariableOp
Conv_1/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_1/conv1d/ExpandDims_1/dimг
Conv_1/conv1d/ExpandDims_1
ExpandDims1Conv_1/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_1/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
Conv_1/conv1d/ExpandDims_1д
Conv_1/conv1dConv2D!Conv_1/conv1d/ExpandDims:output:0#Conv_1/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:џџџџџџџџџА *
paddingVALID*
strides
2
Conv_1/conv1dЈ
Conv_1/conv1d/SqueezeSqueezeConv_1/conv1d:output:0*
T0*,
_output_shapes
:џџџџџџџџџА *
squeeze_dims

§џџџџџџџџ2
Conv_1/conv1d/SqueezeЁ
Conv_1/BiasAdd/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_1/BiasAdd/ReadVariableOpЉ
Conv_1/BiasAddBiasAddConv_1/conv1d/Squeeze:output:0%Conv_1/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
Conv_1/BiasAddr
Conv_1/TanhTanhConv_1/BiasAdd:output:0*
T0*,
_output_shapes
:џџџџџџџџџА 2
Conv_1/Tanh~
max_pooling1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
max_pooling1d/ExpandDims/dimЕ
max_pooling1d/ExpandDims
ExpandDimsConv_1/Tanh:y:0%max_pooling1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:џџџџџџџџџА 2
max_pooling1d/ExpandDimsЩ
max_pooling1d/MaxPoolMaxPool!max_pooling1d/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџX *
ksize
*
paddingVALID*
strides
2
max_pooling1d/MaxPoolІ
max_pooling1d/SqueezeSqueezemax_pooling1d/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџX *
squeeze_dims
2
max_pooling1d/Squeeze
dropout/IdentityIdentitymax_pooling1d/Squeeze:output:0*
T0*+
_output_shapes
:џџџџџџџџџX 2
dropout/Identity
Conv_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
Conv_2/conv1d/ExpandDims/dimО
Conv_2/conv1d/ExpandDims
ExpandDimsdropout/Identity:output:0%Conv_2/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџX 2
Conv_2/conv1d/ExpandDimsЭ
)Conv_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02+
)Conv_2/conv1d/ExpandDims_1/ReadVariableOp
Conv_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_2/conv1d/ExpandDims_1/dimг
Conv_2/conv1d/ExpandDims_1
ExpandDims1Conv_2/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_2/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
Conv_2/conv1d/ExpandDims_1г
Conv_2/conv1dConv2D!Conv_2/conv1d/ExpandDims:output:0#Conv_2/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@*
paddingVALID*
strides
2
Conv_2/conv1dЇ
Conv_2/conv1d/SqueezeSqueezeConv_2/conv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@*
squeeze_dims

§џџџџџџџџ2
Conv_2/conv1d/SqueezeЁ
Conv_2/BiasAdd/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
Conv_2/BiasAdd/ReadVariableOpЈ
Conv_2/BiasAddBiasAddConv_2/conv1d/Squeeze:output:0%Conv_2/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
Conv_2/BiasAddq
Conv_2/TanhTanhConv_2/BiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџL@2
Conv_2/Tanh
max_pooling1d_1/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_1/ExpandDims/dimК
max_pooling1d_1/ExpandDims
ExpandDimsConv_2/Tanh:y:0'max_pooling1d_1/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџL@2
max_pooling1d_1/ExpandDimsЯ
max_pooling1d_1/MaxPoolMaxPool#max_pooling1d_1/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџ&@*
ksize
*
paddingVALID*
strides
2
max_pooling1d_1/MaxPoolЌ
max_pooling1d_1/SqueezeSqueeze max_pooling1d_1/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@*
squeeze_dims
2
max_pooling1d_1/Squeeze
dropout_1/IdentityIdentity max_pooling1d_1/Squeeze:output:0*
T0*+
_output_shapes
:џџџџџџџџџ&@2
dropout_1/Identity
Conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
Conv_3/conv1d/ExpandDims/dimР
Conv_3/conv1d/ExpandDims
ExpandDimsdropout_1/Identity:output:0%Conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ&@2
Conv_3/conv1d/ExpandDimsЭ
)Conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02+
)Conv_3/conv1d/ExpandDims_1/ReadVariableOp
Conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_3/conv1d/ExpandDims_1/dimг
Conv_3/conv1d/ExpandDims_1
ExpandDims1Conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
Conv_3/conv1d/ExpandDims_1г
Conv_3/conv1dConv2D!Conv_3/conv1d/ExpandDims:output:0#Conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  *
paddingVALID*
strides
2
Conv_3/conv1dЇ
Conv_3/conv1d/SqueezeSqueezeConv_3/conv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  *
squeeze_dims

§џџџџџџџџ2
Conv_3/conv1d/SqueezeЁ
Conv_3/BiasAdd/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_3/BiasAdd/ReadVariableOpЈ
Conv_3/BiasAddBiasAddConv_3/conv1d/Squeeze:output:0%Conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
Conv_3/BiasAddq
Conv_3/TanhTanhConv_3/BiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
Conv_3/Tanh
max_pooling1d_2/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_2/ExpandDims/dimК
max_pooling1d_2/ExpandDims
ExpandDimsConv_3/Tanh:y:0'max_pooling1d_2/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  2
max_pooling1d_2/ExpandDimsЯ
max_pooling1d_2/MaxPoolMaxPool#max_pooling1d_2/ExpandDims:output:0*/
_output_shapes
:џџџџџџџџџ *
ksize
*
paddingVALID*
strides
2
max_pooling1d_2/MaxPoolЌ
max_pooling1d_2/SqueezeSqueeze max_pooling1d_2/MaxPool:output:0*
T0*+
_output_shapes
:џџџџџџџџџ *
squeeze_dims
2
max_pooling1d_2/Squeeze
dropout_2/IdentityIdentity max_pooling1d_2/Squeeze:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout_2/Identityo
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
flatten/Const
flatten/ReshapeReshapedropout_2/Identity:output:0flatten/Const:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2
flatten/ReshapeІ
Dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	@*
dtype02
Dense_1/MatMul/ReadVariableOp
Dense_1/MatMulMatMulflatten/Reshape:output:0%Dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
Dense_1/MatMulЄ
Dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02 
Dense_1/BiasAdd/ReadVariableOpЁ
Dense_1/BiasAddBiasAddDense_1/MatMul:product:0&Dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
Dense_1/BiasAddp
Dense_1/TanhTanhDense_1/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ@2
Dense_1/TanhЅ
Dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype02
Dense_2/MatMul/ReadVariableOp
Dense_2/MatMulMatMulDense_1/Tanh:y:0%Dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_2/MatMulЄ
Dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
Dense_2/BiasAdd/ReadVariableOpЁ
Dense_2/BiasAddBiasAddDense_2/MatMul:product:0&Dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_2/BiasAddp
Dense_2/TanhTanhDense_2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_2/TanhЬ
*Dense_3_with_Softmax/MatMul/ReadVariableOpReadVariableOp3dense_3_with_softmax_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*Dense_3_with_Softmax/MatMul/ReadVariableOpМ
Dense_3_with_Softmax/MatMulMatMulDense_2/Tanh:y:02Dense_3_with_Softmax/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_3_with_Softmax/MatMulЫ
+Dense_3_with_Softmax/BiasAdd/ReadVariableOpReadVariableOp4dense_3_with_softmax_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+Dense_3_with_Softmax/BiasAdd/ReadVariableOpе
Dense_3_with_Softmax/BiasAddBiasAdd%Dense_3_with_Softmax/MatMul:product:03Dense_3_with_Softmax/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_3_with_Softmax/BiasAdd 
Dense_3_with_Softmax/SoftmaxSoftmax%Dense_3_with_Softmax/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
Dense_3_with_Softmax/Softmaxй
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulС
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mulй
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulС
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mulй
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulС
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mulЬ
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulФ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulЫ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulФ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulђ
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOp3dense_3_with_softmax_matmul_readvariableop_resource*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulы
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOp4dense_3_with_softmax_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЇ	
IdentityIdentity&Dense_3_with_Softmax/Softmax:softmax:0^Conv_1/BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp*^Conv_1/conv1d/ExpandDims_1/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp*^Conv_2/conv1d/ExpandDims_1/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp*^Conv_3/conv1d/ExpandDims_1/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp^Dense_1/BiasAdd/ReadVariableOp^Dense_1/MatMul/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^Dense_2/BiasAdd/ReadVariableOp^Dense_2/MatMul/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp,^Dense_3_with_Softmax/BiasAdd/ReadVariableOp+^Dense_3_with_Softmax/MatMul/ReadVariableOp<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2>
Conv_1/BiasAdd/ReadVariableOpConv_1/BiasAdd/ReadVariableOp2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2V
)Conv_1/conv1d/ExpandDims_1/ReadVariableOp)Conv_1/conv1d/ExpandDims_1/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2>
Conv_2/BiasAdd/ReadVariableOpConv_2/BiasAdd/ReadVariableOp2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2V
)Conv_2/conv1d/ExpandDims_1/ReadVariableOp)Conv_2/conv1d/ExpandDims_1/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2>
Conv_3/BiasAdd/ReadVariableOpConv_3/BiasAdd/ReadVariableOp2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2V
)Conv_3/conv1d/ExpandDims_1/ReadVariableOp)Conv_3/conv1d/ExpandDims_1/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2@
Dense_1/BiasAdd/ReadVariableOpDense_1/BiasAdd/ReadVariableOp2>
Dense_1/MatMul/ReadVariableOpDense_1/MatMul/ReadVariableOp2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2@
Dense_2/BiasAdd/ReadVariableOpDense_2/BiasAdd/ReadVariableOp2>
Dense_2/MatMul/ReadVariableOpDense_2/MatMul/ReadVariableOp2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2Z
+Dense_3_with_Softmax/BiasAdd/ReadVariableOp+Dense_3_with_Softmax/BiasAdd/ReadVariableOp2X
*Dense_3_with_Softmax/MatMul/ReadVariableOp*Dense_3_with_Softmax/MatMul/ReadVariableOp2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
вЌ
Ь

J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919407

inputs#
conv_1_919142: 
conv_1_919144: #
conv_2_919184: @
conv_2_919186:@#
conv_3_919226:@ 
conv_3_919228: !
dense_1_919271:	@
dense_1_919273:@ 
dense_2_919300:@
dense_2_919302:-
dense_3_with_softmax_919329:)
dense_3_with_softmax_919331:
identityЂConv_1/StatefulPartitionedCallЂ-Conv_1/bias/Regularizer/Square/ReadVariableOpЂ/Conv_1/kernel/Regularizer/Square/ReadVariableOpЂConv_2/StatefulPartitionedCallЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpЂConv_3/StatefulPartitionedCallЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂDense_1/StatefulPartitionedCallЂ.Dense_1/bias/Regularizer/Square/ReadVariableOpЂ0Dense_1/kernel/Regularizer/Square/ReadVariableOpЂDense_2/StatefulPartitionedCallЂ.Dense_2/bias/Regularizer/Square/ReadVariableOpЂ0Dense_2/kernel/Regularizer/Square/ReadVariableOpЂ,Dense_3_with_Softmax/StatefulPartitionedCallЂ;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpЂ=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp
Conv_1/StatefulPartitionedCallStatefulPartitionedCallinputsconv_1_919142conv_1_919144*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџА *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9191412 
Conv_1/StatefulPartitionedCall
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9190702
max_pooling1d/PartitionedCallѕ
dropout/PartitionedCallPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9191532
dropout/PartitionedCallЈ
Conv_2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0conv_2_919184conv_2_919186*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџL@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9191832 
Conv_2/StatefulPartitionedCall
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9190852!
max_pooling1d_1/PartitionedCall§
dropout_1/PartitionedCallPartitionedCall(max_pooling1d_1/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9191952
dropout_1/PartitionedCallЊ
Conv_3/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0conv_3_919226conv_3_919228*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9192252 
Conv_3/StatefulPartitionedCall
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9191002!
max_pooling1d_2/PartitionedCall§
dropout_2/PartitionedCallPartitionedCall(max_pooling1d_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9192372
dropout_2/PartitionedCallю
flatten/PartitionedCallPartitionedCall"dropout_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:џџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9192452
flatten/PartitionedCallЉ
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_919271dense_1_919273*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9192702!
Dense_1/StatefulPartitionedCallБ
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_919300dense_2_919302*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9192992!
Dense_2/StatefulPartitionedCallђ
,Dense_3_with_Softmax/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_with_softmax_919329dense_3_with_softmax_919331*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *Y
fTRR
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_9193282.
,Dense_3_with_Softmax/StatefulPartitionedCallД
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919142*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOpД
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/ConstЖ
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_1/kernel/Regularizer/mul/xИ
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mulЈ
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_919144*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOpІ
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/ConstЎ
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_1/bias/Regularizer/mul/xА
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mulД
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919184*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mulЈ
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_919186*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mulД
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919226*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulЈ
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_919228*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mulД
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919271*
_output_shapes
:	@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOpД
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	@2#
!Dense_1/kernel/Regularizer/Square
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/ConstК
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_1/kernel/Regularizer/mul/xМ
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mulЋ
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_919273*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOpЉ
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/ConstВ
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_1/bias/Regularizer/mul/xД
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mulГ
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919300*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOpГ
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/ConstК
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2"
 Dense_2/kernel/Regularizer/mul/xМ
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mulЋ
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_919302*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOpЉ
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/ConstВ
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2 
Dense_2/bias/Regularizer/mul/xД
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mulк
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919329*
_output_shapes

:*
dtype02?
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOpк
.Dense_3_with_Softmax/kernel/Regularizer/SquareSquareEDense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:20
.Dense_3_with_Softmax/kernel/Regularizer/SquareЏ
-Dense_3_with_Softmax/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2/
-Dense_3_with_Softmax/kernel/Regularizer/Constю
+Dense_3_with_Softmax/kernel/Regularizer/SumSum2Dense_3_with_Softmax/kernel/Regularizer/Square:y:06Dense_3_with_Softmax/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/SumЃ
-Dense_3_with_Softmax/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2/
-Dense_3_with_Softmax/kernel/Regularizer/mul/x№
+Dense_3_with_Softmax/kernel/Regularizer/mulMul6Dense_3_with_Softmax/kernel/Regularizer/mul/x:output:04Dense_3_with_Softmax/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2-
+Dense_3_with_Softmax/kernel/Regularizer/mulв
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_with_softmax_919331*
_output_shapes
:*
dtype02=
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOpа
,Dense_3_with_Softmax/bias/Regularizer/SquareSquareCDense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2.
,Dense_3_with_Softmax/bias/Regularizer/SquareЄ
+Dense_3_with_Softmax/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+Dense_3_with_Softmax/bias/Regularizer/Constц
)Dense_3_with_Softmax/bias/Regularizer/SumSum0Dense_3_with_Softmax/bias/Regularizer/Square:y:04Dense_3_with_Softmax/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/Sum
+Dense_3_with_Softmax/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2-
+Dense_3_with_Softmax/bias/Regularizer/mul/xш
)Dense_3_with_Softmax/bias/Regularizer/mulMul4Dense_3_with_Softmax/bias/Regularizer/mul/x:output:02Dense_3_with_Softmax/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2+
)Dense_3_with_Softmax/bias/Regularizer/mulЫ
IdentityIdentity5Dense_3_with_Softmax/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp-^Dense_3_with_Softmax/StatefulPartitionedCall<^Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp>^Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 2@
Conv_1/StatefulPartitionedCallConv_1/StatefulPartitionedCall2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2@
Conv_2/StatefulPartitionedCallConv_2/StatefulPartitionedCall2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2@
Conv_3/StatefulPartitionedCallConv_3/StatefulPartitionedCall2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2B
Dense_1/StatefulPartitionedCallDense_1/StatefulPartitionedCall2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2B
Dense_2/StatefulPartitionedCallDense_2/StatefulPartitionedCall2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2\
,Dense_3_with_Softmax/StatefulPartitionedCall,Dense_3_with_Softmax/StatefulPartitionedCall2z
;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp;Dense_3_with_Softmax/bias/Regularizer/Square/ReadVariableOp2~
=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp=Dense_3_with_Softmax/kernel/Regularizer/Square/ReadVariableOp:T P
,
_output_shapes
:џџџџџџџџџШ
 
_user_specified_nameinputs
Н
ё
"__inference__traced_restore_921281
file_prefix4
assignvariableop_conv_1_kernel: ,
assignvariableop_1_conv_1_bias: 6
 assignvariableop_2_conv_2_kernel: @,
assignvariableop_3_conv_2_bias:@6
 assignvariableop_4_conv_3_kernel:@ ,
assignvariableop_5_conv_3_bias: 4
!assignvariableop_6_dense_1_kernel:	@-
assignvariableop_7_dense_1_bias:@3
!assignvariableop_8_dense_2_kernel:@-
assignvariableop_9_dense_2_bias:A
/assignvariableop_10_dense_3_with_softmax_kernel:;
-assignvariableop_11_dense_3_with_softmax_bias:'
assignvariableop_12_adam_iter:	 )
assignvariableop_13_adam_beta_1: )
assignvariableop_14_adam_beta_2: (
assignvariableop_15_adam_decay: #
assignvariableop_16_total: #
assignvariableop_17_count: %
assignvariableop_18_total_1: %
assignvariableop_19_count_1: >
(assignvariableop_20_adam_conv_1_kernel_m: 4
&assignvariableop_21_adam_conv_1_bias_m: >
(assignvariableop_22_adam_conv_2_kernel_m: @4
&assignvariableop_23_adam_conv_2_bias_m:@>
(assignvariableop_24_adam_conv_3_kernel_m:@ 4
&assignvariableop_25_adam_conv_3_bias_m: <
)assignvariableop_26_adam_dense_1_kernel_m:	@5
'assignvariableop_27_adam_dense_1_bias_m:@;
)assignvariableop_28_adam_dense_2_kernel_m:@5
'assignvariableop_29_adam_dense_2_bias_m:H
6assignvariableop_30_adam_dense_3_with_softmax_kernel_m:B
4assignvariableop_31_adam_dense_3_with_softmax_bias_m:>
(assignvariableop_32_adam_conv_1_kernel_v: 4
&assignvariableop_33_adam_conv_1_bias_v: >
(assignvariableop_34_adam_conv_2_kernel_v: @4
&assignvariableop_35_adam_conv_2_bias_v:@>
(assignvariableop_36_adam_conv_3_kernel_v:@ 4
&assignvariableop_37_adam_conv_3_bias_v: <
)assignvariableop_38_adam_dense_1_kernel_v:	@5
'assignvariableop_39_adam_dense_1_bias_v:@;
)assignvariableop_40_adam_dense_2_kernel_v:@5
'assignvariableop_41_adam_dense_2_bias_v:H
6assignvariableop_42_adam_dense_3_with_softmax_kernel_v:B
4assignvariableop_43_adam_dense_3_with_softmax_bias_v:
identity_45ЂAssignVariableOpЂAssignVariableOp_1ЂAssignVariableOp_10ЂAssignVariableOp_11ЂAssignVariableOp_12ЂAssignVariableOp_13ЂAssignVariableOp_14ЂAssignVariableOp_15ЂAssignVariableOp_16ЂAssignVariableOp_17ЂAssignVariableOp_18ЂAssignVariableOp_19ЂAssignVariableOp_2ЂAssignVariableOp_20ЂAssignVariableOp_21ЂAssignVariableOp_22ЂAssignVariableOp_23ЂAssignVariableOp_24ЂAssignVariableOp_25ЂAssignVariableOp_26ЂAssignVariableOp_27ЂAssignVariableOp_28ЂAssignVariableOp_29ЂAssignVariableOp_3ЂAssignVariableOp_30ЂAssignVariableOp_31ЂAssignVariableOp_32ЂAssignVariableOp_33ЂAssignVariableOp_34ЂAssignVariableOp_35ЂAssignVariableOp_36ЂAssignVariableOp_37ЂAssignVariableOp_38ЂAssignVariableOp_39ЂAssignVariableOp_4ЂAssignVariableOp_40ЂAssignVariableOp_41ЂAssignVariableOp_42ЂAssignVariableOp_43ЂAssignVariableOp_5ЂAssignVariableOp_6ЂAssignVariableOp_7ЂAssignVariableOp_8ЂAssignVariableOp_9
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*
valueB-B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_namesш
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*m
valuedBb-B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*Ъ
_output_shapesЗ
Д:::::::::::::::::::::::::::::::::::::::::::::*;
dtypes1
/2-	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

Identity
AssignVariableOpAssignVariableOpassignvariableop_conv_1_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1Ѓ
AssignVariableOp_1AssignVariableOpassignvariableop_1_conv_1_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2Ѕ
AssignVariableOp_2AssignVariableOp assignvariableop_2_conv_2_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3Ѓ
AssignVariableOp_3AssignVariableOpassignvariableop_3_conv_2_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4Ѕ
AssignVariableOp_4AssignVariableOp assignvariableop_4_conv_3_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5Ѓ
AssignVariableOp_5AssignVariableOpassignvariableop_5_conv_3_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6І
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_1_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7Є
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_1_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8І
AssignVariableOp_8AssignVariableOp!assignvariableop_8_dense_2_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9Є
AssignVariableOp_9AssignVariableOpassignvariableop_9_dense_2_biasIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10З
AssignVariableOp_10AssignVariableOp/assignvariableop_10_dense_3_with_softmax_kernelIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11Е
AssignVariableOp_11AssignVariableOp-assignvariableop_11_dense_3_with_softmax_biasIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_12Ѕ
AssignVariableOp_12AssignVariableOpassignvariableop_12_adam_iterIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13Ї
AssignVariableOp_13AssignVariableOpassignvariableop_13_adam_beta_1Identity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14Ї
AssignVariableOp_14AssignVariableOpassignvariableop_14_adam_beta_2Identity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15І
AssignVariableOp_15AssignVariableOpassignvariableop_15_adam_decayIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16Ё
AssignVariableOp_16AssignVariableOpassignvariableop_16_totalIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17Ё
AssignVariableOp_17AssignVariableOpassignvariableop_17_countIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18Ѓ
AssignVariableOp_18AssignVariableOpassignvariableop_18_total_1Identity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19Ѓ
AssignVariableOp_19AssignVariableOpassignvariableop_19_count_1Identity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20А
AssignVariableOp_20AssignVariableOp(assignvariableop_20_adam_conv_1_kernel_mIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21Ў
AssignVariableOp_21AssignVariableOp&assignvariableop_21_adam_conv_1_bias_mIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22А
AssignVariableOp_22AssignVariableOp(assignvariableop_22_adam_conv_2_kernel_mIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23Ў
AssignVariableOp_23AssignVariableOp&assignvariableop_23_adam_conv_2_bias_mIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24А
AssignVariableOp_24AssignVariableOp(assignvariableop_24_adam_conv_3_kernel_mIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25Ў
AssignVariableOp_25AssignVariableOp&assignvariableop_25_adam_conv_3_bias_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26Б
AssignVariableOp_26AssignVariableOp)assignvariableop_26_adam_dense_1_kernel_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27Џ
AssignVariableOp_27AssignVariableOp'assignvariableop_27_adam_dense_1_bias_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28Б
AssignVariableOp_28AssignVariableOp)assignvariableop_28_adam_dense_2_kernel_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29Џ
AssignVariableOp_29AssignVariableOp'assignvariableop_29_adam_dense_2_bias_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30О
AssignVariableOp_30AssignVariableOp6assignvariableop_30_adam_dense_3_with_softmax_kernel_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31М
AssignVariableOp_31AssignVariableOp4assignvariableop_31_adam_dense_3_with_softmax_bias_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32А
AssignVariableOp_32AssignVariableOp(assignvariableop_32_adam_conv_1_kernel_vIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33Ў
AssignVariableOp_33AssignVariableOp&assignvariableop_33_adam_conv_1_bias_vIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34А
AssignVariableOp_34AssignVariableOp(assignvariableop_34_adam_conv_2_kernel_vIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35Ў
AssignVariableOp_35AssignVariableOp&assignvariableop_35_adam_conv_2_bias_vIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36А
AssignVariableOp_36AssignVariableOp(assignvariableop_36_adam_conv_3_kernel_vIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37Ў
AssignVariableOp_37AssignVariableOp&assignvariableop_37_adam_conv_3_bias_vIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38Б
AssignVariableOp_38AssignVariableOp)assignvariableop_38_adam_dense_1_kernel_vIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39Џ
AssignVariableOp_39AssignVariableOp'assignvariableop_39_adam_dense_1_bias_vIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40Б
AssignVariableOp_40AssignVariableOp)assignvariableop_40_adam_dense_2_kernel_vIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41Џ
AssignVariableOp_41AssignVariableOp'assignvariableop_41_adam_dense_2_bias_vIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42О
AssignVariableOp_42AssignVariableOp6assignvariableop_42_adam_dense_3_with_softmax_kernel_vIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43М
AssignVariableOp_43AssignVariableOp4assignvariableop_43_adam_dense_3_with_softmax_bias_vIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_439
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOpІ
Identity_44Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_44
Identity_45IdentityIdentity_44:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_45"#
identity_45Identity_45:output:0*m
_input_shapes\
Z: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix

J
.__inference_max_pooling1d_layer_call_fn_919076

inputs
identityн
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9190702
PartitionedCall
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ:e a
=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
Ф
D
(__inference_flatten_layer_call_fn_920714

inputs
identityТ
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:џџџџџџџџџ* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9192452
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
Ю
F
*__inference_dropout_2_layer_call_fn_920687

inputs
identityЧ
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9192372
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs

g
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_919085

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

ExpandDimsБ
MaxPoolMaxPoolExpandDims:output:0*A
_output_shapes/
-:+џџџџџџџџџџџџџџџџџџџџџџџџџџџ*
ksize
*
paddingVALID*
strides
2	
MaxPool
SqueezeSqueezeMaxPool:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ:e a
=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
я
Ѕ
__inference_loss_fn_3_920896D
6conv_2_bias_regularizer_square_readvariableop_resource:@
identityЂ-Conv_2/bias/Regularizer/Square/ReadVariableOpб
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp6conv_2_bias_regularizer_square_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOpІ
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/ConstЎ
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_2/bias/Regularizer/mul/xА
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul
IdentityIdentityConv_2/bias/Regularizer/mul:z:0.^Conv_2/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp

e
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_919070

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

ExpandDimsБ
MaxPoolMaxPoolExpandDims:output:0*A
_output_shapes/
-:+џџџџџџџџџџџџџџџџџџџџџџџџџџџ*
ksize
*
paddingVALID*
strides
2	
MaxPool
SqueezeSqueezeMaxPool:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ:e a
=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
а$
ѓ
B__inference_Conv_3_layer_call_and_return_conditional_losses_920682

inputsA
+conv1d_expanddims_1_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identityЂBiasAdd/ReadVariableOpЂ-Conv_3/bias/Regularizer/Square/ReadVariableOpЂ/Conv_3/kernel/Regularizer/Square/ReadVariableOpЂ"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
§џџџџџџџџ2
conv1d/ExpandDims/dim
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:џџџџџџџџџ&@2
conv1d/ExpandDimsИ
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dimЗ
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
conv1d/ExpandDims_1З
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:џџџџџџџџџ  *
paddingVALID*
strides
2
conv1d
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  *
squeeze_dims

§џџџџџџџџ2
conv1d/Squeeze
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:џџџџџџџџџ  2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:џџџџџџџџџ  2
Tanhв
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOpД
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/ConstЖ
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_3/kernel/Regularizer/mul/xИ
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mulК
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOpІ
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/ConstЎ
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2
Conv_3/bias/Regularizer/mul/xА
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:џџџџџџџџџ  2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:џџџџџџџџџ&@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:џџџџџџџџџ&@
 
_user_specified_nameinputs
к
c
*__inference_dropout_2_layer_call_fn_920692

inputs
identityЂStatefulPartitionedCallп
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9194902
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
к

Л
$__inference_signature_wrapper_920102
convolutional_inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identityЂStatefulPartitionedCallт
StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8 **
f%R#
!__inference__wrapped_model_9190612
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:џџџџџџџџџШ: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
,
_output_shapes
:џџџџџџџџџШ
.
_user_specified_nameConvolutional_inputs
Ь
d
E__inference_dropout_2_layer_call_and_return_conditional_losses_919490

inputs
identityc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *ЂМ?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/ShapeИ
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:џџџџџџџџџ *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЭЬL=2
dropout/GreaterEqual/yТ
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:џџџџџџџџџ 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:џџџџџџџџџ 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:џџџџџџџџџ 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ :S O
+
_output_shapes
:џџџџџџџџџ 
 
_user_specified_nameinputs
Т
Б
__inference_loss_fn_2_920885N
8conv_2_kernel_regularizer_square_readvariableop_resource: @
identityЂ/Conv_2/kernel/Regularizer/Square/ReadVariableOpп
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp8conv_2_kernel_regularizer_square_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOpД
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/ConstЖ
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o:2!
Conv_2/kernel/Regularizer/mul/xИ
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul
IdentityIdentity!Conv_2/kernel/Regularizer/mul:z:00^Conv_2/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp

g
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_919100

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

ExpandDimsБ
MaxPoolMaxPoolExpandDims:output:0*A
_output_shapes/
-:+џџџџџџџџџџџџџџџџџџџџџџџџџџџ*
ksize
*
paddingVALID*
strides
2	
MaxPool
SqueezeSqueezeMaxPool:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ:e a
=
_output_shapes+
):'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
ж
a
(__inference_dropout_layer_call_fn_920540

inputs
identityЂStatefulPartitionedCallн
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџX * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9195562
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:џџџџџџџџџX 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџX 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:џџџџџџџџџX 
 
_user_specified_nameinputs"ЬL
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*ж
serving_defaultТ
Z
Convolutional_inputsB
&serving_default_Convolutional_inputs:0џџџџџџџџџШH
Dense_3_with_Softmax0
StatefulPartitionedCall:0џџџџџџџџџtensorflow/serving/predict:ов
Џ~
layer-0
layer_with_weights-0
layer-1
layer-2
layer-3
layer_with_weights-1
layer-4
layer-5
layer-6
layer_with_weights-2
layer-7
	layer-8

layer-9
layer-10
layer_with_weights-3
layer-11
layer_with_weights-4
layer-12
layer_with_weights-5
layer-13
	optimizer
	variables
regularization_losses
trainable_variables
	keras_api

signatures
Т__call__
+У&call_and_return_all_conditional_losses
Ф_default_save_signature"z
_tf_keras_networkяy{"name": "FCN_classifier", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "class_name": "Functional", "config": {"name": "FCN_classifier", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "Convolutional_inputs"}, "name": "Convolutional_inputs", "inbound_nodes": []}, {"class_name": "Conv1D", "config": {"name": "Conv_1", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [25]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_1", "inbound_nodes": [[["Convolutional_inputs", 0, 0, {}]]]}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d", "inbound_nodes": [[["Conv_1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["max_pooling1d", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "Conv_2", "trainable": true, "dtype": "float32", "filters": 64, "kernel_size": {"class_name": "__tuple__", "items": [13]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_2", "inbound_nodes": [[["dropout", 0, 0, {}]]]}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_1", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_1", "inbound_nodes": [[["Conv_2", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout_1", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "name": "dropout_1", "inbound_nodes": [[["max_pooling1d_1", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "Conv_3", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_3", "inbound_nodes": [[["dropout_1", 0, 0, {}]]]}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_2", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_2", "inbound_nodes": [[["Conv_3", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "name": "dropout_2", "inbound_nodes": [[["max_pooling1d_2", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["dropout_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "Dense_1", "trainable": true, "dtype": "float32", "units": 64, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_1", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "Dense_2", "trainable": true, "dtype": "float32", "units": 16, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_2", "inbound_nodes": [[["Dense_1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "Dense_3_with_Softmax", "trainable": true, "dtype": "float32", "units": 8, "activation": "softmax", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_3_with_Softmax", "inbound_nodes": [[["Dense_2", 0, 0, {}]]]}], "input_layers": [["Convolutional_inputs", 0, 0]], "output_layers": [["Dense_3_with_Softmax", 0, 0]]}, "shared_object_id": 38, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 200, 8]}, "is_graph_network": true, "save_spec": {"class_name": "TypeSpec", "type_spec": "tf.TensorSpec", "serialized": [{"class_name": "TensorShape", "items": [null, 200, 8]}, "float32", "Convolutional_inputs"]}, "keras_version": "2.5.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "FCN_classifier", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "Convolutional_inputs"}, "name": "Convolutional_inputs", "inbound_nodes": [], "shared_object_id": 0}, {"class_name": "Conv1D", "config": {"name": "Conv_1", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [25]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 1}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 2}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 3}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 4}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_1", "inbound_nodes": [[["Convolutional_inputs", 0, 0, {}]]], "shared_object_id": 5}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d", "inbound_nodes": [[["Conv_1", 0, 0, {}]]], "shared_object_id": 6}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["max_pooling1d", 0, 0, {}]]], "shared_object_id": 7}, {"class_name": "Conv1D", "config": {"name": "Conv_2", "trainable": true, "dtype": "float32", "filters": 64, "kernel_size": {"class_name": "__tuple__", "items": [13]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 8}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 9}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 10}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 11}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_2", "inbound_nodes": [[["dropout", 0, 0, {}]]], "shared_object_id": 12}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_1", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_1", "inbound_nodes": [[["Conv_2", 0, 0, {}]]], "shared_object_id": 13}, {"class_name": "Dropout", "config": {"name": "dropout_1", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "name": "dropout_1", "inbound_nodes": [[["max_pooling1d_1", 0, 0, {}]]], "shared_object_id": 14}, {"class_name": "Conv1D", "config": {"name": "Conv_3", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 15}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 16}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 17}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 18}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_3", "inbound_nodes": [[["dropout_1", 0, 0, {}]]], "shared_object_id": 19}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_2", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_2", "inbound_nodes": [[["Conv_3", 0, 0, {}]]], "shared_object_id": 20}, {"class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "name": "dropout_2", "inbound_nodes": [[["max_pooling1d_2", 0, 0, {}]]], "shared_object_id": 21}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["dropout_2", 0, 0, {}]]], "shared_object_id": 22}, {"class_name": "Dense", "config": {"name": "Dense_1", "trainable": true, "dtype": "float32", "units": 64, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 23}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 24}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 25}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 26}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_1", "inbound_nodes": [[["flatten", 0, 0, {}]]], "shared_object_id": 27}, {"class_name": "Dense", "config": {"name": "Dense_2", "trainable": true, "dtype": "float32", "units": 16, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 28}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 29}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 30}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 31}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_2", "inbound_nodes": [[["Dense_1", 0, 0, {}]]], "shared_object_id": 32}, {"class_name": "Dense", "config": {"name": "Dense_3_with_Softmax", "trainable": true, "dtype": "float32", "units": 8, "activation": "softmax", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 33}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 34}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 35}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 36}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_3_with_Softmax", "inbound_nodes": [[["Dense_2", 0, 0, {}]]], "shared_object_id": 37}], "input_layers": [["Convolutional_inputs", 0, 0]], "output_layers": [["Dense_3_with_Softmax", 0, 0]]}}, "training_config": {"loss": "sparse_categorical_crossentropy", "metrics": [[{"class_name": "MeanMetricWrapper", "config": {"name": "accuracy", "dtype": "float32", "fn": "sparse_categorical_accuracy"}, "shared_object_id": 40}]], "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": {"class_name": "CosineDecay", "config": {"initial_learning_rate": 0.001, "decay_steps": 50000, "alpha": 0.05, "name": null}, "shared_object_id": 41}, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
"
_tf_keras_input_layerъ{"class_name": "InputLayer", "name": "Convolutional_inputs", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "Convolutional_inputs"}}


kernel
bias
	variables
regularization_losses
trainable_variables
	keras_api
Х__call__
+Ц&call_and_return_all_conditional_losses"ѓ

_tf_keras_layerй
{"name": "Conv_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Conv1D", "config": {"name": "Conv_1", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [25]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 1}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 2}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 3}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 4}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["Convolutional_inputs", 0, 0, {}]]], "shared_object_id": 5, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 8}}, "shared_object_id": 42}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 200, 8]}}
б
	variables
regularization_losses
trainable_variables
	keras_api
Ч__call__
+Ш&call_and_return_all_conditional_losses"Р
_tf_keras_layerІ{"name": "max_pooling1d", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "MaxPooling1D", "config": {"name": "max_pooling1d", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "inbound_nodes": [[["Conv_1", 0, 0, {}]]], "shared_object_id": 6, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}, "shared_object_id": 43}}
­
	variables
 regularization_losses
!trainable_variables
"	keras_api
Щ__call__
+Ъ&call_and_return_all_conditional_losses"
_tf_keras_layer{"name": "dropout", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "inbound_nodes": [[["max_pooling1d", 0, 0, {}]]], "shared_object_id": 7}


#kernel
$bias
%	variables
&regularization_losses
'trainable_variables
(	keras_api
Ы__call__
+Ь&call_and_return_all_conditional_losses"ъ

_tf_keras_layerа
{"name": "Conv_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Conv1D", "config": {"name": "Conv_2", "trainable": true, "dtype": "float32", "filters": 64, "kernel_size": {"class_name": "__tuple__", "items": [13]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 8}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 9}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 10}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 11}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["dropout", 0, 0, {}]]], "shared_object_id": 12, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 32}}, "shared_object_id": 44}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 88, 32]}}
ж
)	variables
*regularization_losses
+trainable_variables
,	keras_api
Э__call__
+Ю&call_and_return_all_conditional_losses"Х
_tf_keras_layerЋ{"name": "max_pooling1d_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_1", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "inbound_nodes": [[["Conv_2", 0, 0, {}]]], "shared_object_id": 13, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}, "shared_object_id": 45}}
Д
-	variables
.regularization_losses
/trainable_variables
0	keras_api
Я__call__
+а&call_and_return_all_conditional_losses"Ѓ
_tf_keras_layer{"name": "dropout_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dropout", "config": {"name": "dropout_1", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "inbound_nodes": [[["max_pooling1d_1", 0, 0, {}]]], "shared_object_id": 14}


1kernel
2bias
3	variables
4regularization_losses
5trainable_variables
6	keras_api
б__call__
+в&call_and_return_all_conditional_losses"э

_tf_keras_layerг
{"name": "Conv_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Conv1D", "config": {"name": "Conv_3", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 15}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 16}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 17}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 18}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["dropout_1", 0, 0, {}]]], "shared_object_id": 19, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 64}}, "shared_object_id": 46}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 38, 64]}}
ж
7	variables
8regularization_losses
9trainable_variables
:	keras_api
г__call__
+д&call_and_return_all_conditional_losses"Х
_tf_keras_layerЋ{"name": "max_pooling1d_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_2", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "inbound_nodes": [[["Conv_3", 0, 0, {}]]], "shared_object_id": 20, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}, "shared_object_id": 47}}
Д
;	variables
<regularization_losses
=trainable_variables
>	keras_api
е__call__
+ж&call_and_return_all_conditional_losses"Ѓ
_tf_keras_layer{"name": "dropout_2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.05, "noise_shape": null, "seed": null}, "inbound_nodes": [[["max_pooling1d_2", 0, 0, {}]]], "shared_object_id": 21}
Т
?	variables
@regularization_losses
Atrainable_variables
B	keras_api
з__call__
+и&call_and_return_all_conditional_losses"Б
_tf_keras_layer{"name": "flatten", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "inbound_nodes": [[["dropout_2", 0, 0, {}]]], "shared_object_id": 22, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}, "shared_object_id": 48}}
Ђ


Ckernel
Dbias
E	variables
Fregularization_losses
Gtrainable_variables
H	keras_api
й__call__
+к&call_and_return_all_conditional_losses"ћ
_tf_keras_layerс{"name": "Dense_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dense", "config": {"name": "Dense_1", "trainable": true, "dtype": "float32", "units": 64, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 23}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 24}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 25}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 26}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["flatten", 0, 0, {}]]], "shared_object_id": 27, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 512}}, "shared_object_id": 49}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 512]}}
 


Ikernel
Jbias
K	variables
Lregularization_losses
Mtrainable_variables
N	keras_api
л__call__
+м&call_and_return_all_conditional_losses"љ
_tf_keras_layerп{"name": "Dense_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dense", "config": {"name": "Dense_2", "trainable": true, "dtype": "float32", "units": 16, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 28}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 29}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 30}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 31}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["Dense_1", 0, 0, {}]]], "shared_object_id": 32, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 64}}, "shared_object_id": 50}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 64]}}
М


Okernel
Pbias
Q	variables
Rregularization_losses
Strainable_variables
T	keras_api
н__call__
+о&call_and_return_all_conditional_losses"	
_tf_keras_layerћ{"name": "Dense_3_with_Softmax", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dense", "config": {"name": "Dense_3_with_Softmax", "trainable": true, "dtype": "float32", "units": 8, "activation": "softmax", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 33}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 34}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 35}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 36}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["Dense_2", 0, 0, {}]]], "shared_object_id": 37, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 16}}, "shared_object_id": 51}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 16]}}
А
Uiter

Vbeta_1

Wbeta_2
	XdecaymЊmЋ#mЌ$m­1mЎ2mЏCmАDmБImВJmГOmДPmЕvЖvЗ#vИ$vЙ1vК2vЛCvМDvНIvОJvПOvРPvС"
	optimizer
v
0
1
#2
$3
14
25
C6
D7
I8
J9
O10
P11"
trackable_list_wrapper

п0
р1
с2
т3
у4
ф5
х6
ц7
ч8
ш9
щ10
ъ11"
trackable_list_wrapper
v
0
1
#2
$3
14
25
C6
D7
I8
J9
O10
P11"
trackable_list_wrapper
Ю
Ylayer_regularization_losses
	variables
Zmetrics

[layers
regularization_losses
trainable_variables
\non_trainable_variables
]layer_metrics
Т__call__
Ф_default_save_signature
+У&call_and_return_all_conditional_losses
'У"call_and_return_conditional_losses"
_generic_user_object
-
ыserving_default"
signature_map
#:! 2Conv_1/kernel
: 2Conv_1/bias
.
0
1"
trackable_list_wrapper
0
п0
р1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
А
^layer_regularization_losses
	variables
_metrics

`layers
regularization_losses
trainable_variables
anon_trainable_variables
blayer_metrics
Х__call__
+Ц&call_and_return_all_conditional_losses
'Ц"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
А
clayer_regularization_losses
	variables
dmetrics

elayers
regularization_losses
trainable_variables
fnon_trainable_variables
glayer_metrics
Ч__call__
+Ш&call_and_return_all_conditional_losses
'Ш"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
А
hlayer_regularization_losses
	variables
imetrics

jlayers
 regularization_losses
!trainable_variables
knon_trainable_variables
llayer_metrics
Щ__call__
+Ъ&call_and_return_all_conditional_losses
'Ъ"call_and_return_conditional_losses"
_generic_user_object
#:! @2Conv_2/kernel
:@2Conv_2/bias
.
#0
$1"
trackable_list_wrapper
0
с0
т1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
А
mlayer_regularization_losses
%	variables
nmetrics

olayers
&regularization_losses
'trainable_variables
pnon_trainable_variables
qlayer_metrics
Ы__call__
+Ь&call_and_return_all_conditional_losses
'Ь"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
А
rlayer_regularization_losses
)	variables
smetrics

tlayers
*regularization_losses
+trainable_variables
unon_trainable_variables
vlayer_metrics
Э__call__
+Ю&call_and_return_all_conditional_losses
'Ю"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
А
wlayer_regularization_losses
-	variables
xmetrics

ylayers
.regularization_losses
/trainable_variables
znon_trainable_variables
{layer_metrics
Я__call__
+а&call_and_return_all_conditional_losses
'а"call_and_return_conditional_losses"
_generic_user_object
#:!@ 2Conv_3/kernel
: 2Conv_3/bias
.
10
21"
trackable_list_wrapper
0
у0
ф1"
trackable_list_wrapper
.
10
21"
trackable_list_wrapper
Б
|layer_regularization_losses
3	variables
}metrics

~layers
4regularization_losses
5trainable_variables
non_trainable_variables
layer_metrics
б__call__
+в&call_and_return_all_conditional_losses
'в"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
Е
 layer_regularization_losses
7	variables
metrics
layers
8regularization_losses
9trainable_variables
non_trainable_variables
layer_metrics
г__call__
+д&call_and_return_all_conditional_losses
'д"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
Е
 layer_regularization_losses
;	variables
metrics
layers
<regularization_losses
=trainable_variables
non_trainable_variables
layer_metrics
е__call__
+ж&call_and_return_all_conditional_losses
'ж"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
Е
 layer_regularization_losses
?	variables
metrics
layers
@regularization_losses
Atrainable_variables
non_trainable_variables
layer_metrics
з__call__
+и&call_and_return_all_conditional_losses
'и"call_and_return_conditional_losses"
_generic_user_object
!:	@2Dense_1/kernel
:@2Dense_1/bias
.
C0
D1"
trackable_list_wrapper
0
х0
ц1"
trackable_list_wrapper
.
C0
D1"
trackable_list_wrapper
Е
 layer_regularization_losses
E	variables
metrics
layers
Fregularization_losses
Gtrainable_variables
non_trainable_variables
layer_metrics
й__call__
+к&call_and_return_all_conditional_losses
'к"call_and_return_conditional_losses"
_generic_user_object
 :@2Dense_2/kernel
:2Dense_2/bias
.
I0
J1"
trackable_list_wrapper
0
ч0
ш1"
trackable_list_wrapper
.
I0
J1"
trackable_list_wrapper
Е
 layer_regularization_losses
K	variables
metrics
layers
Lregularization_losses
Mtrainable_variables
non_trainable_variables
layer_metrics
л__call__
+м&call_and_return_all_conditional_losses
'м"call_and_return_conditional_losses"
_generic_user_object
-:+2Dense_3_with_Softmax/kernel
':%2Dense_3_with_Softmax/bias
.
O0
P1"
trackable_list_wrapper
0
щ0
ъ1"
trackable_list_wrapper
.
O0
P1"
trackable_list_wrapper
Е
 layer_regularization_losses
Q	variables
metrics
layers
Rregularization_losses
Strainable_variables
non_trainable_variables
layer_metrics
н__call__
+о&call_and_return_all_conditional_losses
'о"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
 "
trackable_list_wrapper
0
0
 1"
trackable_list_wrapper

0
1
2
3
4
5
6
7
	8

9
10
11
12
13"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
п0
р1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
с0
т1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
у0
ф1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
х0
ц1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
ч0
ш1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
щ0
ъ1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
и

Ёtotal

Ђcount
Ѓ	variables
Є	keras_api"
_tf_keras_metric{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}, "shared_object_id": 52}
Ѓ

Ѕtotal

Іcount
Ї
_fn_kwargs
Ј	variables
Љ	keras_api"з
_tf_keras_metricМ{"class_name": "MeanMetricWrapper", "name": "accuracy", "dtype": "float32", "config": {"name": "accuracy", "dtype": "float32", "fn": "sparse_categorical_accuracy"}, "shared_object_id": 40}
:  (2total
:  (2count
0
Ё0
Ђ1"
trackable_list_wrapper
.
Ѓ	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
Ѕ0
І1"
trackable_list_wrapper
.
Ј	variables"
_generic_user_object
(:& 2Adam/Conv_1/kernel/m
: 2Adam/Conv_1/bias/m
(:& @2Adam/Conv_2/kernel/m
:@2Adam/Conv_2/bias/m
(:&@ 2Adam/Conv_3/kernel/m
: 2Adam/Conv_3/bias/m
&:$	@2Adam/Dense_1/kernel/m
:@2Adam/Dense_1/bias/m
%:#@2Adam/Dense_2/kernel/m
:2Adam/Dense_2/bias/m
2:02"Adam/Dense_3_with_Softmax/kernel/m
,:*2 Adam/Dense_3_with_Softmax/bias/m
(:& 2Adam/Conv_1/kernel/v
: 2Adam/Conv_1/bias/v
(:& @2Adam/Conv_2/kernel/v
:@2Adam/Conv_2/bias/v
(:&@ 2Adam/Conv_3/kernel/v
: 2Adam/Conv_3/bias/v
&:$	@2Adam/Dense_1/kernel/v
:@2Adam/Dense_1/bias/v
%:#@2Adam/Dense_2/kernel/v
:2Adam/Dense_2/bias/v
2:02"Adam/Dense_3_with_Softmax/kernel/v
,:*2 Adam/Dense_3_with_Softmax/bias/v
2
/__inference_FCN_classifier_layer_call_fn_919434
/__inference_FCN_classifier_layer_call_fn_920131
/__inference_FCN_classifier_layer_call_fn_920160
/__inference_FCN_classifier_layer_call_fn_919769Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
і2ѓ
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_920310
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_920481
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919882
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919995Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ё2ю
!__inference__wrapped_model_919061Ш
В
FullArgSpec
args 
varargsjargs
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *8Ђ5
30
Convolutional_inputsџџџџџџџџџШ
б2Ю
'__inference_Conv_1_layer_call_fn_920502Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ь2щ
B__inference_Conv_1_layer_call_and_return_conditional_losses_920530Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
2
.__inference_max_pooling1d_layer_call_fn_919076г
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *3Ђ0
.+'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Є2Ё
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_919070г
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *3Ђ0
.+'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
2
(__inference_dropout_layer_call_fn_920535
(__inference_dropout_layer_call_fn_920540Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
Ф2С
C__inference_dropout_layer_call_and_return_conditional_losses_920545
C__inference_dropout_layer_call_and_return_conditional_losses_920557Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
б2Ю
'__inference_Conv_2_layer_call_fn_920578Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ь2щ
B__inference_Conv_2_layer_call_and_return_conditional_losses_920606Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
2
0__inference_max_pooling1d_1_layer_call_fn_919091г
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *3Ђ0
.+'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
І2Ѓ
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_919085г
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *3Ђ0
.+'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
2
*__inference_dropout_1_layer_call_fn_920611
*__inference_dropout_1_layer_call_fn_920616Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
Ш2Х
E__inference_dropout_1_layer_call_and_return_conditional_losses_920621
E__inference_dropout_1_layer_call_and_return_conditional_losses_920633Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
б2Ю
'__inference_Conv_3_layer_call_fn_920654Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ь2щ
B__inference_Conv_3_layer_call_and_return_conditional_losses_920682Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
2
0__inference_max_pooling1d_2_layer_call_fn_919106г
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *3Ђ0
.+'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
І2Ѓ
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_919100г
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *3Ђ0
.+'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
2
*__inference_dropout_2_layer_call_fn_920687
*__inference_dropout_2_layer_call_fn_920692Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
Ш2Х
E__inference_dropout_2_layer_call_and_return_conditional_losses_920697
E__inference_dropout_2_layer_call_and_return_conditional_losses_920709Д
ЋВЇ
FullArgSpec)
args!
jself
jinputs

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
в2Я
(__inference_flatten_layer_call_fn_920714Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
э2ъ
C__inference_flatten_layer_call_and_return_conditional_losses_920720Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
в2Я
(__inference_Dense_1_layer_call_fn_920741Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
э2ъ
C__inference_Dense_1_layer_call_and_return_conditional_losses_920764Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
в2Я
(__inference_Dense_2_layer_call_fn_920785Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
э2ъ
C__inference_Dense_2_layer_call_and_return_conditional_losses_920808Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
п2м
5__inference_Dense_3_with_Softmax_layer_call_fn_920829Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
њ2ї
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_920852Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
Г2А
__inference_loss_fn_0_920863
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_1_920874
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_2_920885
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_3_920896
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_4_920907
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_5_920918
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_6_920929
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_7_920940
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_8_920951
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Г2А
__inference_loss_fn_9_920962
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Д2Б
__inference_loss_fn_10_920973
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
Д2Б
__inference_loss_fn_11_920984
В
FullArgSpec
args 
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *Ђ 
иBе
$__inference_signature_wrapper_920102Convolutional_inputs"
В
FullArgSpec
args 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 Ќ
B__inference_Conv_1_layer_call_and_return_conditional_losses_920530f4Ђ1
*Ђ'
%"
inputsџџџџџџџџџШ
Њ "*Ђ'
 
0џџџџџџџџџА 
 
'__inference_Conv_1_layer_call_fn_920502Y4Ђ1
*Ђ'
%"
inputsџџџџџџџџџШ
Њ "џџџџџџџџџА Њ
B__inference_Conv_2_layer_call_and_return_conditional_losses_920606d#$3Ђ0
)Ђ&
$!
inputsџџџџџџџџџX 
Њ ")Ђ&

0џџџџџџџџџL@
 
'__inference_Conv_2_layer_call_fn_920578W#$3Ђ0
)Ђ&
$!
inputsџџџџџџџџџX 
Њ "џџџџџџџџџL@Њ
B__inference_Conv_3_layer_call_and_return_conditional_losses_920682d123Ђ0
)Ђ&
$!
inputsџџџџџџџџџ&@
Њ ")Ђ&

0џџџџџџџџџ  
 
'__inference_Conv_3_layer_call_fn_920654W123Ђ0
)Ђ&
$!
inputsџџџџџџџџџ&@
Њ "џџџџџџџџџ  Є
C__inference_Dense_1_layer_call_and_return_conditional_losses_920764]CD0Ђ-
&Ђ#
!
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ@
 |
(__inference_Dense_1_layer_call_fn_920741PCD0Ђ-
&Ђ#
!
inputsџџџџџџџџџ
Њ "џџџџџџџџџ@Ѓ
C__inference_Dense_2_layer_call_and_return_conditional_losses_920808\IJ/Ђ,
%Ђ"
 
inputsџџџџџџџџџ@
Њ "%Ђ"

0џџџџџџџџџ
 {
(__inference_Dense_2_layer_call_fn_920785OIJ/Ђ,
%Ђ"
 
inputsџџџџџџџџџ@
Њ "џџџџџџџџџА
P__inference_Dense_3_with_Softmax_layer_call_and_return_conditional_losses_920852\OP/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ
 
5__inference_Dense_3_with_Softmax_layer_call_fn_920829OOP/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџа
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919882#$12CDIJOPJЂG
@Ђ=
30
Convolutional_inputsџџџџџџџџџШ
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 а
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_919995#$12CDIJOPJЂG
@Ђ=
30
Convolutional_inputsџџџџџџџџџШ
p

 
Њ "%Ђ"

0џџџџџџџџџ
 С
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_920310s#$12CDIJOP<Ђ9
2Ђ/
%"
inputsџџџџџџџџџШ
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 С
J__inference_FCN_classifier_layer_call_and_return_conditional_losses_920481s#$12CDIJOP<Ђ9
2Ђ/
%"
inputsџџџџџџџџџШ
p

 
Њ "%Ђ"

0џџџџџџџџџ
 Ї
/__inference_FCN_classifier_layer_call_fn_919434t#$12CDIJOPJЂG
@Ђ=
30
Convolutional_inputsџџџџџџџџџШ
p 

 
Њ "џџџџџџџџџЇ
/__inference_FCN_classifier_layer_call_fn_919769t#$12CDIJOPJЂG
@Ђ=
30
Convolutional_inputsџџџџџџџџџШ
p

 
Њ "џџџџџџџџџ
/__inference_FCN_classifier_layer_call_fn_920131f#$12CDIJOP<Ђ9
2Ђ/
%"
inputsџџџџџџџџџШ
p 

 
Њ "џџџџџџџџџ
/__inference_FCN_classifier_layer_call_fn_920160f#$12CDIJOP<Ђ9
2Ђ/
%"
inputsџџџџџџџџџШ
p

 
Њ "џџџџџџџџџХ
!__inference__wrapped_model_919061#$12CDIJOPBЂ?
8Ђ5
30
Convolutional_inputsџџџџџџџџџШ
Њ "KЊH
F
Dense_3_with_Softmax.+
Dense_3_with_Softmaxџџџџџџџџџ­
E__inference_dropout_1_layer_call_and_return_conditional_losses_920621d7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ&@
p 
Њ ")Ђ&

0џџџџџџџџџ&@
 ­
E__inference_dropout_1_layer_call_and_return_conditional_losses_920633d7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ&@
p
Њ ")Ђ&

0џџџџџџџџџ&@
 
*__inference_dropout_1_layer_call_fn_920611W7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ&@
p 
Њ "џџџџџџџџџ&@
*__inference_dropout_1_layer_call_fn_920616W7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ&@
p
Њ "џџџџџџџџџ&@­
E__inference_dropout_2_layer_call_and_return_conditional_losses_920697d7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ 
p 
Њ ")Ђ&

0џџџџџџџџџ 
 ­
E__inference_dropout_2_layer_call_and_return_conditional_losses_920709d7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ 
p
Њ ")Ђ&

0џџџџџџџџџ 
 
*__inference_dropout_2_layer_call_fn_920687W7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ 
p 
Њ "џџџџџџџџџ 
*__inference_dropout_2_layer_call_fn_920692W7Ђ4
-Ђ*
$!
inputsџџџџџџџџџ 
p
Њ "џџџџџџџџџ Ћ
C__inference_dropout_layer_call_and_return_conditional_losses_920545d7Ђ4
-Ђ*
$!
inputsџџџџџџџџџX 
p 
Њ ")Ђ&

0џџџџџџџџџX 
 Ћ
C__inference_dropout_layer_call_and_return_conditional_losses_920557d7Ђ4
-Ђ*
$!
inputsџџџџџџџџџX 
p
Њ ")Ђ&

0џџџџџџџџџX 
 
(__inference_dropout_layer_call_fn_920535W7Ђ4
-Ђ*
$!
inputsџџџџџџџџџX 
p 
Њ "џџџџџџџџџX 
(__inference_dropout_layer_call_fn_920540W7Ђ4
-Ђ*
$!
inputsџџџџџџџџџX 
p
Њ "џџџџџџџџџX Є
C__inference_flatten_layer_call_and_return_conditional_losses_920720]3Ђ0
)Ђ&
$!
inputsџџџџџџџџџ 
Њ "&Ђ#

0џџџџџџџџџ
 |
(__inference_flatten_layer_call_fn_920714P3Ђ0
)Ђ&
$!
inputsџџџџџџџџџ 
Њ "џџџџџџџџџ;
__inference_loss_fn_0_920863Ђ

Ђ 
Њ " <
__inference_loss_fn_10_920973OЂ

Ђ 
Њ " <
__inference_loss_fn_11_920984PЂ

Ђ 
Њ " ;
__inference_loss_fn_1_920874Ђ

Ђ 
Њ " ;
__inference_loss_fn_2_920885#Ђ

Ђ 
Њ " ;
__inference_loss_fn_3_920896$Ђ

Ђ 
Њ " ;
__inference_loss_fn_4_9209071Ђ

Ђ 
Њ " ;
__inference_loss_fn_5_9209182Ђ

Ђ 
Њ " ;
__inference_loss_fn_6_920929CЂ

Ђ 
Њ " ;
__inference_loss_fn_7_920940DЂ

Ђ 
Њ " ;
__inference_loss_fn_8_920951IЂ

Ђ 
Њ " ;
__inference_loss_fn_9_920962JЂ

Ђ 
Њ " д
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_919085EЂB
;Ђ8
63
inputs'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Њ ";Ђ8
1.
0'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 Ћ
0__inference_max_pooling1d_1_layer_call_fn_919091wEЂB
;Ђ8
63
inputs'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Њ ".+'џџџџџџџџџџџџџџџџџџџџџџџџџџџд
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_919100EЂB
;Ђ8
63
inputs'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Њ ";Ђ8
1.
0'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 Ћ
0__inference_max_pooling1d_2_layer_call_fn_919106wEЂB
;Ђ8
63
inputs'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Њ ".+'џџџџџџџџџџџџџџџџџџџџџџџџџџџв
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_919070EЂB
;Ђ8
63
inputs'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Њ ";Ђ8
1.
0'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
 Љ
.__inference_max_pooling1d_layer_call_fn_919076wEЂB
;Ђ8
63
inputs'џџџџџџџџџџџџџџџџџџџџџџџџџџџ
Њ ".+'џџџџџџџџџџџџџџџџџџџџџџџџџџџр
$__inference_signature_wrapper_920102З#$12CDIJOPZЂW
Ђ 
PЊM
K
Convolutional_inputs30
Convolutional_inputsџџџџџџџџџШ"KЊH
F
Dense_3_with_Softmax.+
Dense_3_with_Softmaxџџџџџџџџџ