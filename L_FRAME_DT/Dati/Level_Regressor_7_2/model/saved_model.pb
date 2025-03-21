��
��
B
AssignVariableOp
resource
value"dtype"
dtypetype�
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
�
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
�
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
delete_old_dirsbool(�
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
dtypetype�
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
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
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
N
Squeeze

input"T
output"T"	
Ttype"
squeeze_dims	list(int)
 (
�
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
executor_typestring �
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
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.5.02v2.5.0-rc3-213-ga4dfb8d1a718��
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
shape:	�@*
shared_nameDense_1/kernel
r
"Dense_1/kernel/Read/ReadVariableOpReadVariableOpDense_1/kernel*
_output_shapes
:	�@*
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
x
Dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*
shared_nameDense_3/kernel
q
"Dense_3/kernel/Read/ReadVariableOpReadVariableOpDense_3/kernel*
_output_shapes

:*
dtype0
p
Dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameDense_3/bias
i
 Dense_3/bias/Read/ReadVariableOpReadVariableOpDense_3/bias*
_output_shapes
:*
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
�
Adam/Conv_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/Conv_1/kernel/m
�
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
�
Adam/Conv_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*%
shared_nameAdam/Conv_2/kernel/m
�
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
�
Adam/Conv_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:@ *%
shared_nameAdam/Conv_3/kernel/m
�
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
�
Adam/Dense_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�@*&
shared_nameAdam/Dense_1/kernel/m
�
)Adam/Dense_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Dense_1/kernel/m*
_output_shapes
:	�@*
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
�
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
�
Adam/Dense_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/Dense_3/kernel/m

)Adam/Dense_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/Dense_3/kernel/m*
_output_shapes

:*
dtype0
~
Adam/Dense_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/Dense_3/bias/m
w
'Adam/Dense_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/Dense_3/bias/m*
_output_shapes
:*
dtype0
�
Adam/Conv_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/Conv_1/kernel/v
�
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
�
Adam/Conv_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*%
shared_nameAdam/Conv_2/kernel/v
�
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
�
Adam/Conv_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:@ *%
shared_nameAdam/Conv_3/kernel/v
�
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
�
Adam/Dense_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�@*&
shared_nameAdam/Dense_1/kernel/v
�
)Adam/Dense_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Dense_1/kernel/v*
_output_shapes
:	�@*
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
�
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
�
Adam/Dense_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/Dense_3/kernel/v

)Adam/Dense_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/Dense_3/kernel/v*
_output_shapes

:*
dtype0
~
Adam/Dense_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/Dense_3/bias/v
w
'Adam/Dense_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/Dense_3/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
�O
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�N
value�NB�N B�N
�
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
�
Uiter

Vbeta_1

Wbeta_2
	Xdecaym�m�#m�$m�1m�2m�Cm�Dm�Im�Jm�Om�Pm�v�v�#v�$v�1v�2v�Cv�Dv�Iv�Jv�Ov�Pv�
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
�
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
�
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
�
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
�
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
�
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
�
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
�
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
�
|layer_regularization_losses
3	variables
}metrics

~layers
4regularization_losses
5trainable_variables
non_trainable_variables
�layer_metrics
 
 
 
�
 �layer_regularization_losses
7	variables
�metrics
�layers
8regularization_losses
9trainable_variables
�non_trainable_variables
�layer_metrics
 
 
 
�
 �layer_regularization_losses
;	variables
�metrics
�layers
<regularization_losses
=trainable_variables
�non_trainable_variables
�layer_metrics
 
 
 
�
 �layer_regularization_losses
?	variables
�metrics
�layers
@regularization_losses
Atrainable_variables
�non_trainable_variables
�layer_metrics
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
�
 �layer_regularization_losses
E	variables
�metrics
�layers
Fregularization_losses
Gtrainable_variables
�non_trainable_variables
�layer_metrics
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
�
 �layer_regularization_losses
K	variables
�metrics
�layers
Lregularization_losses
Mtrainable_variables
�non_trainable_variables
�layer_metrics
ZX
VARIABLE_VALUEDense_3/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEDense_3/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE

O0
P1
 

O0
P1
�
 �layer_regularization_losses
Q	variables
�metrics
�layers
Rregularization_losses
Strainable_variables
�non_trainable_variables
�layer_metrics
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
�0
�1
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

�total

�count
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
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
}{
VARIABLE_VALUEAdam/Dense_3/kernel/mRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/Dense_3/bias/mPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
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
}{
VARIABLE_VALUEAdam/Dense_3/kernel/vRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/Dense_3/bias/vPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
�
$serving_default_Convolutional_inputsPlaceholder*,
_output_shapes
:����������*
dtype0*!
shape:����������
�
StatefulPartitionedCallStatefulPartitionedCall$serving_default_Convolutional_inputsConv_1/kernelConv_1/biasConv_2/kernelConv_2/biasConv_3/kernelConv_3/biasDense_1/kernelDense_1/biasDense_2/kernelDense_2/biasDense_3/kernelDense_3/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *-
f(R&
$__inference_signature_wrapper_984962
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename!Conv_1/kernel/Read/ReadVariableOpConv_1/bias/Read/ReadVariableOp!Conv_2/kernel/Read/ReadVariableOpConv_2/bias/Read/ReadVariableOp!Conv_3/kernel/Read/ReadVariableOpConv_3/bias/Read/ReadVariableOp"Dense_1/kernel/Read/ReadVariableOp Dense_1/bias/Read/ReadVariableOp"Dense_2/kernel/Read/ReadVariableOp Dense_2/bias/Read/ReadVariableOp"Dense_3/kernel/Read/ReadVariableOp Dense_3/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp(Adam/Conv_1/kernel/m/Read/ReadVariableOp&Adam/Conv_1/bias/m/Read/ReadVariableOp(Adam/Conv_2/kernel/m/Read/ReadVariableOp&Adam/Conv_2/bias/m/Read/ReadVariableOp(Adam/Conv_3/kernel/m/Read/ReadVariableOp&Adam/Conv_3/bias/m/Read/ReadVariableOp)Adam/Dense_1/kernel/m/Read/ReadVariableOp'Adam/Dense_1/bias/m/Read/ReadVariableOp)Adam/Dense_2/kernel/m/Read/ReadVariableOp'Adam/Dense_2/bias/m/Read/ReadVariableOp)Adam/Dense_3/kernel/m/Read/ReadVariableOp'Adam/Dense_3/bias/m/Read/ReadVariableOp(Adam/Conv_1/kernel/v/Read/ReadVariableOp&Adam/Conv_1/bias/v/Read/ReadVariableOp(Adam/Conv_2/kernel/v/Read/ReadVariableOp&Adam/Conv_2/bias/v/Read/ReadVariableOp(Adam/Conv_3/kernel/v/Read/ReadVariableOp&Adam/Conv_3/bias/v/Read/ReadVariableOp)Adam/Dense_1/kernel/v/Read/ReadVariableOp'Adam/Dense_1/bias/v/Read/ReadVariableOp)Adam/Dense_2/kernel/v/Read/ReadVariableOp'Adam/Dense_2/bias/v/Read/ReadVariableOp)Adam/Dense_3/kernel/v/Read/ReadVariableOp'Adam/Dense_3/bias/v/Read/ReadVariableOpConst*9
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
GPU 2J 8� *(
f#R!
__inference__traced_save_985926
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameConv_1/kernelConv_1/biasConv_2/kernelConv_2/biasConv_3/kernelConv_3/biasDense_1/kernelDense_1/biasDense_2/kernelDense_2/biasDense_3/kernelDense_3/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decaytotalcounttotal_1count_1Adam/Conv_1/kernel/mAdam/Conv_1/bias/mAdam/Conv_2/kernel/mAdam/Conv_2/bias/mAdam/Conv_3/kernel/mAdam/Conv_3/bias/mAdam/Dense_1/kernel/mAdam/Dense_1/bias/mAdam/Dense_2/kernel/mAdam/Dense_2/bias/mAdam/Dense_3/kernel/mAdam/Dense_3/bias/mAdam/Conv_1/kernel/vAdam/Conv_1/bias/vAdam/Conv_2/kernel/vAdam/Conv_2/bias/vAdam/Conv_3/kernel/vAdam/Conv_3/bias/vAdam/Dense_1/kernel/vAdam/Dense_1/bias/vAdam/Dense_2/kernel/vAdam/Dense_2/bias/vAdam/Dense_3/kernel/vAdam/Dense_3/bias/v*8
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
GPU 2J 8� *+
f&R$
"__inference__traced_restore_986068�
�
c
*__inference_dropout_1_layer_call_fn_985450

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9844312
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:���������&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������&@22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�
g
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_984018

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim�

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+���������������������������2

ExpandDims�
MaxPoolMaxPoolExpandDims:output:0*A
_output_shapes/
-:+���������������������������*
ksize
*
paddingVALID*
strides
2	
MaxPool�
SqueezeSqueezeMaxPool:output:0*
T0*=
_output_shapes+
):'���������������������������*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'���������������������������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'���������������������������:e a
=
_output_shapes+
):'���������������������������
 
_user_specified_nameinputs
�
c
*__inference_dropout_2_layer_call_fn_985526

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9843982
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:��������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
��
�	
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984315

inputs#
conv_1_984075: 
conv_1_984077: #
conv_2_984117: @
conv_2_984119:@#
conv_3_984159:@ 
conv_3_984161: !
dense_1_984204:	�@
dense_1_984206:@ 
dense_2_984233:@
dense_2_984235: 
dense_3_984249:
dense_3_984251:
identity��Conv_1/StatefulPartitionedCall�-Conv_1/bias/Regularizer/Square/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�Conv_2/StatefulPartitionedCall�-Conv_2/bias/Regularizer/Square/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�Conv_3/StatefulPartitionedCall�-Conv_3/bias/Regularizer/Square/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�Dense_1/StatefulPartitionedCall�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�Dense_2/StatefulPartitionedCall�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�Dense_3/StatefulPartitionedCall�
Conv_1/StatefulPartitionedCallStatefulPartitionedCallinputsconv_1_984075conv_1_984077*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:���������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9840742 
Conv_1/StatefulPartitionedCall�
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9840032
max_pooling1d/PartitionedCall�
dropout/PartitionedCallPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9840862
dropout/PartitionedCall�
Conv_2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0conv_2_984117conv_2_984119*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������L@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9841162 
Conv_2/StatefulPartitionedCall�
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9840182!
max_pooling1d_1/PartitionedCall�
dropout_1/PartitionedCallPartitionedCall(max_pooling1d_1/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9841282
dropout_1/PartitionedCall�
Conv_3/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0conv_3_984159conv_3_984161*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9841582 
Conv_3/StatefulPartitionedCall�
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9840332!
max_pooling1d_2/PartitionedCall�
dropout_2/PartitionedCallPartitionedCall(max_pooling1d_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9841702
dropout_2/PartitionedCall�
flatten/PartitionedCallPartitionedCall"dropout_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9841782
flatten/PartitionedCall�
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_984204dense_1_984206*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9842032!
Dense_1/StatefulPartitionedCall�
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_984233dense_2_984235*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9842322!
Dense_2/StatefulPartitionedCall�
Dense_3/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_984249dense_3_984251*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_3_layer_call_and_return_conditional_losses_9842482!
Dense_3/StatefulPartitionedCall�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984075*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984077*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984117*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984119*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984159*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984161*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984204*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984206*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984233*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984235*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentity(Dense_3/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp ^Dense_3/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2@
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
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2B
Dense_3/StatefulPartitionedCallDense_3/StatefulPartitionedCall:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_4_985716N
8conv_3_kernel_regularizer_square_readvariableop_resource:@ 
identity��/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp8conv_3_kernel_regularizer_square_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
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
�
d
E__inference_dropout_2_layer_call_and_return_conditional_losses_984398

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:��������� 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:��������� *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:��������� 2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:��������� 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:��������� 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:��������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�	
�
C__inference_Dense_3_layer_call_and_return_conditional_losses_984248

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdd�
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
L
0__inference_max_pooling1d_2_layer_call_fn_984039

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *=
_output_shapes+
):'���������������������������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9840332
PartitionedCall�
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'���������������������������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'���������������������������:e a
=
_output_shapes+
):'���������������������������
 
_user_specified_nameinputs
�
_
C__inference_flatten_layer_call_and_return_conditional_losses_984178

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"����   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:����������2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
��
�
!__inference__wrapped_model_983994
convolutional_inputsV
@fcn_regressor_conv_1_conv1d_expanddims_1_readvariableop_resource: B
4fcn_regressor_conv_1_biasadd_readvariableop_resource: V
@fcn_regressor_conv_2_conv1d_expanddims_1_readvariableop_resource: @B
4fcn_regressor_conv_2_biasadd_readvariableop_resource:@V
@fcn_regressor_conv_3_conv1d_expanddims_1_readvariableop_resource:@ B
4fcn_regressor_conv_3_biasadd_readvariableop_resource: G
4fcn_regressor_dense_1_matmul_readvariableop_resource:	�@C
5fcn_regressor_dense_1_biasadd_readvariableop_resource:@F
4fcn_regressor_dense_2_matmul_readvariableop_resource:@C
5fcn_regressor_dense_2_biasadd_readvariableop_resource:F
4fcn_regressor_dense_3_matmul_readvariableop_resource:C
5fcn_regressor_dense_3_biasadd_readvariableop_resource:
identity��+FCN_regressor/Conv_1/BiasAdd/ReadVariableOp�7FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOp�+FCN_regressor/Conv_2/BiasAdd/ReadVariableOp�7FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOp�+FCN_regressor/Conv_3/BiasAdd/ReadVariableOp�7FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOp�,FCN_regressor/Dense_1/BiasAdd/ReadVariableOp�+FCN_regressor/Dense_1/MatMul/ReadVariableOp�,FCN_regressor/Dense_2/BiasAdd/ReadVariableOp�+FCN_regressor/Dense_2/MatMul/ReadVariableOp�,FCN_regressor/Dense_3/BiasAdd/ReadVariableOp�+FCN_regressor/Dense_3/MatMul/ReadVariableOp�
*FCN_regressor/Conv_1/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2,
*FCN_regressor/Conv_1/conv1d/ExpandDims/dim�
&FCN_regressor/Conv_1/conv1d/ExpandDims
ExpandDimsconvolutional_inputs3FCN_regressor/Conv_1/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:����������2(
&FCN_regressor/Conv_1/conv1d/ExpandDims�
7FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp@fcn_regressor_conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype029
7FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOp�
,FCN_regressor/Conv_1/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2.
,FCN_regressor/Conv_1/conv1d/ExpandDims_1/dim�
(FCN_regressor/Conv_1/conv1d/ExpandDims_1
ExpandDims?FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOp:value:05FCN_regressor/Conv_1/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2*
(FCN_regressor/Conv_1/conv1d/ExpandDims_1�
FCN_regressor/Conv_1/conv1dConv2D/FCN_regressor/Conv_1/conv1d/ExpandDims:output:01FCN_regressor/Conv_1/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:���������� *
paddingVALID*
strides
2
FCN_regressor/Conv_1/conv1d�
#FCN_regressor/Conv_1/conv1d/SqueezeSqueeze$FCN_regressor/Conv_1/conv1d:output:0*
T0*,
_output_shapes
:���������� *
squeeze_dims

���������2%
#FCN_regressor/Conv_1/conv1d/Squeeze�
+FCN_regressor/Conv_1/BiasAdd/ReadVariableOpReadVariableOp4fcn_regressor_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02-
+FCN_regressor/Conv_1/BiasAdd/ReadVariableOp�
FCN_regressor/Conv_1/BiasAddBiasAdd,FCN_regressor/Conv_1/conv1d/Squeeze:output:03FCN_regressor/Conv_1/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:���������� 2
FCN_regressor/Conv_1/BiasAdd�
FCN_regressor/Conv_1/TanhTanh%FCN_regressor/Conv_1/BiasAdd:output:0*
T0*,
_output_shapes
:���������� 2
FCN_regressor/Conv_1/Tanh�
*FCN_regressor/max_pooling1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2,
*FCN_regressor/max_pooling1d/ExpandDims/dim�
&FCN_regressor/max_pooling1d/ExpandDims
ExpandDimsFCN_regressor/Conv_1/Tanh:y:03FCN_regressor/max_pooling1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:���������� 2(
&FCN_regressor/max_pooling1d/ExpandDims�
#FCN_regressor/max_pooling1d/MaxPoolMaxPool/FCN_regressor/max_pooling1d/ExpandDims:output:0*/
_output_shapes
:���������X *
ksize
*
paddingVALID*
strides
2%
#FCN_regressor/max_pooling1d/MaxPool�
#FCN_regressor/max_pooling1d/SqueezeSqueeze,FCN_regressor/max_pooling1d/MaxPool:output:0*
T0*+
_output_shapes
:���������X *
squeeze_dims
2%
#FCN_regressor/max_pooling1d/Squeeze�
FCN_regressor/dropout/IdentityIdentity,FCN_regressor/max_pooling1d/Squeeze:output:0*
T0*+
_output_shapes
:���������X 2 
FCN_regressor/dropout/Identity�
*FCN_regressor/Conv_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2,
*FCN_regressor/Conv_2/conv1d/ExpandDims/dim�
&FCN_regressor/Conv_2/conv1d/ExpandDims
ExpandDims'FCN_regressor/dropout/Identity:output:03FCN_regressor/Conv_2/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������X 2(
&FCN_regressor/Conv_2/conv1d/ExpandDims�
7FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp@fcn_regressor_conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype029
7FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOp�
,FCN_regressor/Conv_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2.
,FCN_regressor/Conv_2/conv1d/ExpandDims_1/dim�
(FCN_regressor/Conv_2/conv1d/ExpandDims_1
ExpandDims?FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOp:value:05FCN_regressor/Conv_2/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2*
(FCN_regressor/Conv_2/conv1d/ExpandDims_1�
FCN_regressor/Conv_2/conv1dConv2D/FCN_regressor/Conv_2/conv1d/ExpandDims:output:01FCN_regressor/Conv_2/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������L@*
paddingVALID*
strides
2
FCN_regressor/Conv_2/conv1d�
#FCN_regressor/Conv_2/conv1d/SqueezeSqueeze$FCN_regressor/Conv_2/conv1d:output:0*
T0*+
_output_shapes
:���������L@*
squeeze_dims

���������2%
#FCN_regressor/Conv_2/conv1d/Squeeze�
+FCN_regressor/Conv_2/BiasAdd/ReadVariableOpReadVariableOp4fcn_regressor_conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02-
+FCN_regressor/Conv_2/BiasAdd/ReadVariableOp�
FCN_regressor/Conv_2/BiasAddBiasAdd,FCN_regressor/Conv_2/conv1d/Squeeze:output:03FCN_regressor/Conv_2/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������L@2
FCN_regressor/Conv_2/BiasAdd�
FCN_regressor/Conv_2/TanhTanh%FCN_regressor/Conv_2/BiasAdd:output:0*
T0*+
_output_shapes
:���������L@2
FCN_regressor/Conv_2/Tanh�
,FCN_regressor/max_pooling1d_1/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2.
,FCN_regressor/max_pooling1d_1/ExpandDims/dim�
(FCN_regressor/max_pooling1d_1/ExpandDims
ExpandDimsFCN_regressor/Conv_2/Tanh:y:05FCN_regressor/max_pooling1d_1/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������L@2*
(FCN_regressor/max_pooling1d_1/ExpandDims�
%FCN_regressor/max_pooling1d_1/MaxPoolMaxPool1FCN_regressor/max_pooling1d_1/ExpandDims:output:0*/
_output_shapes
:���������&@*
ksize
*
paddingVALID*
strides
2'
%FCN_regressor/max_pooling1d_1/MaxPool�
%FCN_regressor/max_pooling1d_1/SqueezeSqueeze.FCN_regressor/max_pooling1d_1/MaxPool:output:0*
T0*+
_output_shapes
:���������&@*
squeeze_dims
2'
%FCN_regressor/max_pooling1d_1/Squeeze�
 FCN_regressor/dropout_1/IdentityIdentity.FCN_regressor/max_pooling1d_1/Squeeze:output:0*
T0*+
_output_shapes
:���������&@2"
 FCN_regressor/dropout_1/Identity�
*FCN_regressor/Conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2,
*FCN_regressor/Conv_3/conv1d/ExpandDims/dim�
&FCN_regressor/Conv_3/conv1d/ExpandDims
ExpandDims)FCN_regressor/dropout_1/Identity:output:03FCN_regressor/Conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������&@2(
&FCN_regressor/Conv_3/conv1d/ExpandDims�
7FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp@fcn_regressor_conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype029
7FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOp�
,FCN_regressor/Conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2.
,FCN_regressor/Conv_3/conv1d/ExpandDims_1/dim�
(FCN_regressor/Conv_3/conv1d/ExpandDims_1
ExpandDims?FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:05FCN_regressor/Conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2*
(FCN_regressor/Conv_3/conv1d/ExpandDims_1�
FCN_regressor/Conv_3/conv1dConv2D/FCN_regressor/Conv_3/conv1d/ExpandDims:output:01FCN_regressor/Conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������  *
paddingVALID*
strides
2
FCN_regressor/Conv_3/conv1d�
#FCN_regressor/Conv_3/conv1d/SqueezeSqueeze$FCN_regressor/Conv_3/conv1d:output:0*
T0*+
_output_shapes
:���������  *
squeeze_dims

���������2%
#FCN_regressor/Conv_3/conv1d/Squeeze�
+FCN_regressor/Conv_3/BiasAdd/ReadVariableOpReadVariableOp4fcn_regressor_conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02-
+FCN_regressor/Conv_3/BiasAdd/ReadVariableOp�
FCN_regressor/Conv_3/BiasAddBiasAdd,FCN_regressor/Conv_3/conv1d/Squeeze:output:03FCN_regressor/Conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������  2
FCN_regressor/Conv_3/BiasAdd�
FCN_regressor/Conv_3/TanhTanh%FCN_regressor/Conv_3/BiasAdd:output:0*
T0*+
_output_shapes
:���������  2
FCN_regressor/Conv_3/Tanh�
,FCN_regressor/max_pooling1d_2/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2.
,FCN_regressor/max_pooling1d_2/ExpandDims/dim�
(FCN_regressor/max_pooling1d_2/ExpandDims
ExpandDimsFCN_regressor/Conv_3/Tanh:y:05FCN_regressor/max_pooling1d_2/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������  2*
(FCN_regressor/max_pooling1d_2/ExpandDims�
%FCN_regressor/max_pooling1d_2/MaxPoolMaxPool1FCN_regressor/max_pooling1d_2/ExpandDims:output:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2'
%FCN_regressor/max_pooling1d_2/MaxPool�
%FCN_regressor/max_pooling1d_2/SqueezeSqueeze.FCN_regressor/max_pooling1d_2/MaxPool:output:0*
T0*+
_output_shapes
:��������� *
squeeze_dims
2'
%FCN_regressor/max_pooling1d_2/Squeeze�
 FCN_regressor/dropout_2/IdentityIdentity.FCN_regressor/max_pooling1d_2/Squeeze:output:0*
T0*+
_output_shapes
:��������� 2"
 FCN_regressor/dropout_2/Identity�
FCN_regressor/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����   2
FCN_regressor/flatten/Const�
FCN_regressor/flatten/ReshapeReshape)FCN_regressor/dropout_2/Identity:output:0$FCN_regressor/flatten/Const:output:0*
T0*(
_output_shapes
:����������2
FCN_regressor/flatten/Reshape�
+FCN_regressor/Dense_1/MatMul/ReadVariableOpReadVariableOp4fcn_regressor_dense_1_matmul_readvariableop_resource*
_output_shapes
:	�@*
dtype02-
+FCN_regressor/Dense_1/MatMul/ReadVariableOp�
FCN_regressor/Dense_1/MatMulMatMul&FCN_regressor/flatten/Reshape:output:03FCN_regressor/Dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
FCN_regressor/Dense_1/MatMul�
,FCN_regressor/Dense_1/BiasAdd/ReadVariableOpReadVariableOp5fcn_regressor_dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02.
,FCN_regressor/Dense_1/BiasAdd/ReadVariableOp�
FCN_regressor/Dense_1/BiasAddBiasAdd&FCN_regressor/Dense_1/MatMul:product:04FCN_regressor/Dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
FCN_regressor/Dense_1/BiasAdd�
FCN_regressor/Dense_1/TanhTanh&FCN_regressor/Dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:���������@2
FCN_regressor/Dense_1/Tanh�
+FCN_regressor/Dense_2/MatMul/ReadVariableOpReadVariableOp4fcn_regressor_dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype02-
+FCN_regressor/Dense_2/MatMul/ReadVariableOp�
FCN_regressor/Dense_2/MatMulMatMulFCN_regressor/Dense_1/Tanh:y:03FCN_regressor/Dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
FCN_regressor/Dense_2/MatMul�
,FCN_regressor/Dense_2/BiasAdd/ReadVariableOpReadVariableOp5fcn_regressor_dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,FCN_regressor/Dense_2/BiasAdd/ReadVariableOp�
FCN_regressor/Dense_2/BiasAddBiasAdd&FCN_regressor/Dense_2/MatMul:product:04FCN_regressor/Dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
FCN_regressor/Dense_2/BiasAdd�
FCN_regressor/Dense_2/TanhTanh&FCN_regressor/Dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
FCN_regressor/Dense_2/Tanh�
+FCN_regressor/Dense_3/MatMul/ReadVariableOpReadVariableOp4fcn_regressor_dense_3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02-
+FCN_regressor/Dense_3/MatMul/ReadVariableOp�
FCN_regressor/Dense_3/MatMulMatMulFCN_regressor/Dense_2/Tanh:y:03FCN_regressor/Dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
FCN_regressor/Dense_3/MatMul�
,FCN_regressor/Dense_3/BiasAdd/ReadVariableOpReadVariableOp5fcn_regressor_dense_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,FCN_regressor/Dense_3/BiasAdd/ReadVariableOp�
FCN_regressor/Dense_3/BiasAddBiasAdd&FCN_regressor/Dense_3/MatMul:product:04FCN_regressor/Dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
FCN_regressor/Dense_3/BiasAdd�
IdentityIdentity&FCN_regressor/Dense_3/BiasAdd:output:0,^FCN_regressor/Conv_1/BiasAdd/ReadVariableOp8^FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOp,^FCN_regressor/Conv_2/BiasAdd/ReadVariableOp8^FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOp,^FCN_regressor/Conv_3/BiasAdd/ReadVariableOp8^FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOp-^FCN_regressor/Dense_1/BiasAdd/ReadVariableOp,^FCN_regressor/Dense_1/MatMul/ReadVariableOp-^FCN_regressor/Dense_2/BiasAdd/ReadVariableOp,^FCN_regressor/Dense_2/MatMul/ReadVariableOp-^FCN_regressor/Dense_3/BiasAdd/ReadVariableOp,^FCN_regressor/Dense_3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2Z
+FCN_regressor/Conv_1/BiasAdd/ReadVariableOp+FCN_regressor/Conv_1/BiasAdd/ReadVariableOp2r
7FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOp7FCN_regressor/Conv_1/conv1d/ExpandDims_1/ReadVariableOp2Z
+FCN_regressor/Conv_2/BiasAdd/ReadVariableOp+FCN_regressor/Conv_2/BiasAdd/ReadVariableOp2r
7FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOp7FCN_regressor/Conv_2/conv1d/ExpandDims_1/ReadVariableOp2Z
+FCN_regressor/Conv_3/BiasAdd/ReadVariableOp+FCN_regressor/Conv_3/BiasAdd/ReadVariableOp2r
7FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOp7FCN_regressor/Conv_3/conv1d/ExpandDims_1/ReadVariableOp2\
,FCN_regressor/Dense_1/BiasAdd/ReadVariableOp,FCN_regressor/Dense_1/BiasAdd/ReadVariableOp2Z
+FCN_regressor/Dense_1/MatMul/ReadVariableOp+FCN_regressor/Dense_1/MatMul/ReadVariableOp2\
,FCN_regressor/Dense_2/BiasAdd/ReadVariableOp,FCN_regressor/Dense_2/BiasAdd/ReadVariableOp2Z
+FCN_regressor/Dense_2/MatMul/ReadVariableOp+FCN_regressor/Dense_2/MatMul/ReadVariableOp2\
,FCN_regressor/Dense_3/BiasAdd/ReadVariableOp,FCN_regressor/Dense_3/BiasAdd/ReadVariableOp2Z
+FCN_regressor/Dense_3/MatMul/ReadVariableOp+FCN_regressor/Dense_3/MatMul/ReadVariableOp:b ^
,
_output_shapes
:����������
.
_user_specified_nameConvolutional_inputs
�
g
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_984033

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim�

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+���������������������������2

ExpandDims�
MaxPoolMaxPoolExpandDims:output:0*A
_output_shapes/
-:+���������������������������*
ksize
*
paddingVALID*
strides
2	
MaxPool�
SqueezeSqueezeMaxPool:output:0*
T0*=
_output_shapes+
):'���������������������������*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'���������������������������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'���������������������������:e a
=
_output_shapes+
):'���������������������������
 
_user_specified_nameinputs
�
�
C__inference_Dense_2_layer_call_and_return_conditional_losses_985642

inputs0
matmul_readvariableop_resource:@-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Tanh�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
a
(__inference_dropout_layer_call_fn_985374

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9844642
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:���������X 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������X 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
��
�

I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984867
convolutional_inputs#
conv_1_984769: 
conv_1_984771: #
conv_2_984776: @
conv_2_984778:@#
conv_3_984783:@ 
conv_3_984785: !
dense_1_984791:	�@
dense_1_984793:@ 
dense_2_984796:@
dense_2_984798: 
dense_3_984801:
dense_3_984803:
identity��Conv_1/StatefulPartitionedCall�-Conv_1/bias/Regularizer/Square/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�Conv_2/StatefulPartitionedCall�-Conv_2/bias/Regularizer/Square/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�Conv_3/StatefulPartitionedCall�-Conv_3/bias/Regularizer/Square/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�Dense_1/StatefulPartitionedCall�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�Dense_2/StatefulPartitionedCall�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�Dense_3/StatefulPartitionedCall�dropout/StatefulPartitionedCall�!dropout_1/StatefulPartitionedCall�!dropout_2/StatefulPartitionedCall�
Conv_1/StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsconv_1_984769conv_1_984771*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:���������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9840742 
Conv_1/StatefulPartitionedCall�
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9840032
max_pooling1d/PartitionedCall�
dropout/StatefulPartitionedCallStatefulPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9844642!
dropout/StatefulPartitionedCall�
Conv_2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0conv_2_984776conv_2_984778*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������L@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9841162 
Conv_2/StatefulPartitionedCall�
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9840182!
max_pooling1d_1/PartitionedCall�
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_1/PartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9844312#
!dropout_1/StatefulPartitionedCall�
Conv_3/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0conv_3_984783conv_3_984785*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9841582 
Conv_3/StatefulPartitionedCall�
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9840332!
max_pooling1d_2/PartitionedCall�
!dropout_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_2/PartitionedCall:output:0"^dropout_1/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9843982#
!dropout_2/StatefulPartitionedCall�
flatten/PartitionedCallPartitionedCall*dropout_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9841782
flatten/PartitionedCall�
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_984791dense_1_984793*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9842032!
Dense_1/StatefulPartitionedCall�
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_984796dense_2_984798*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9842322!
Dense_2/StatefulPartitionedCall�
Dense_3/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_984801dense_3_984803*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_3_layer_call_and_return_conditional_losses_9842482!
Dense_3/StatefulPartitionedCall�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984769*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984771*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984776*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984778*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984783*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984785*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984791*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984793*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984796*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984798*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentity(Dense_3/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp ^Dense_3/StatefulPartitionedCall ^dropout/StatefulPartitionedCall"^dropout_1/StatefulPartitionedCall"^dropout_2/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2@
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
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2B
Dense_3/StatefulPartitionedCallDense_3/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall2F
!dropout_2/StatefulPartitionedCall!dropout_2/StatefulPartitionedCall:b ^
,
_output_shapes
:����������
.
_user_specified_nameConvolutional_inputs
�	
�
C__inference_Dense_3_layer_call_and_return_conditional_losses_985661

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdd�
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
a
C__inference_dropout_layer_call_and_return_conditional_losses_984086

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:���������X 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:���������X 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������X :S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�
�
'__inference_Conv_1_layer_call_fn_985336

inputs
unknown: 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:���������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9840742
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:���������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_8_985760K
9dense_2_kernel_regularizer_square_readvariableop_resource:@
identity��0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_2_kernel_regularizer_square_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
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
�
_
C__inference_flatten_layer_call_and_return_conditional_losses_985554

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"����   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:����������2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
.__inference_FCN_regressor_layer_call_fn_985020

inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	�@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_9846092
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�$
�
B__inference_Conv_1_layer_call_and_return_conditional_losses_984074

inputsA
+conv1d_expanddims_1_readvariableop_resource: -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�-Conv_1/bias/Regularizer/Square/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
conv1d/ExpandDims/dim�
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:����������2
conv1d/ExpandDims�
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
conv1d/ExpandDims_1/dim�
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
conv1d/ExpandDims_1�
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:���������� *
paddingVALID*
strides
2
conv1d�
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:���������� *
squeeze_dims

���������2
conv1d/Squeeze�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:���������� 2	
BiasAdd]
TanhTanhBiasAdd:output:0*
T0*,
_output_shapes
:���������� 2
Tanh�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:���������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
D
(__inference_dropout_layer_call_fn_985369

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9840862
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:���������X 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������X :S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�Z
�
__inference__traced_save_985926
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
'savev2_dense_2_bias_read_readvariableop-
)savev2_dense_3_kernel_read_readvariableop+
'savev2_dense_3_bias_read_readvariableop(
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
.savev2_adam_dense_2_bias_m_read_readvariableop4
0savev2_adam_dense_3_kernel_m_read_readvariableop2
.savev2_adam_dense_3_bias_m_read_readvariableop3
/savev2_adam_conv_1_kernel_v_read_readvariableop1
-savev2_adam_conv_1_bias_v_read_readvariableop3
/savev2_adam_conv_2_kernel_v_read_readvariableop1
-savev2_adam_conv_2_bias_v_read_readvariableop3
/savev2_adam_conv_3_kernel_v_read_readvariableop1
-savev2_adam_conv_3_bias_v_read_readvariableop4
0savev2_adam_dense_1_kernel_v_read_readvariableop2
.savev2_adam_dense_1_bias_v_read_readvariableop4
0savev2_adam_dense_2_kernel_v_read_readvariableop2
.savev2_adam_dense_2_bias_v_read_readvariableop4
0savev2_adam_dense_3_kernel_v_read_readvariableop2
.savev2_adam_dense_3_bias_v_read_readvariableop
savev2_const

identity_1��MergeV2Checkpoints�
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
Const_1�
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
ShardedFilename/shard�
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename�
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*�
value�B�-B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*m
valuedBb-B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices�
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0(savev2_conv_1_kernel_read_readvariableop&savev2_conv_1_bias_read_readvariableop(savev2_conv_2_kernel_read_readvariableop&savev2_conv_2_bias_read_readvariableop(savev2_conv_3_kernel_read_readvariableop&savev2_conv_3_bias_read_readvariableop)savev2_dense_1_kernel_read_readvariableop'savev2_dense_1_bias_read_readvariableop)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop)savev2_dense_3_kernel_read_readvariableop'savev2_dense_3_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop/savev2_adam_conv_1_kernel_m_read_readvariableop-savev2_adam_conv_1_bias_m_read_readvariableop/savev2_adam_conv_2_kernel_m_read_readvariableop-savev2_adam_conv_2_bias_m_read_readvariableop/savev2_adam_conv_3_kernel_m_read_readvariableop-savev2_adam_conv_3_bias_m_read_readvariableop0savev2_adam_dense_1_kernel_m_read_readvariableop.savev2_adam_dense_1_bias_m_read_readvariableop0savev2_adam_dense_2_kernel_m_read_readvariableop.savev2_adam_dense_2_bias_m_read_readvariableop0savev2_adam_dense_3_kernel_m_read_readvariableop.savev2_adam_dense_3_bias_m_read_readvariableop/savev2_adam_conv_1_kernel_v_read_readvariableop-savev2_adam_conv_1_bias_v_read_readvariableop/savev2_adam_conv_2_kernel_v_read_readvariableop-savev2_adam_conv_2_bias_v_read_readvariableop/savev2_adam_conv_3_kernel_v_read_readvariableop-savev2_adam_conv_3_bias_v_read_readvariableop0savev2_adam_dense_1_kernel_v_read_readvariableop.savev2_adam_dense_1_bias_v_read_readvariableop0savev2_adam_dense_2_kernel_v_read_readvariableop.savev2_adam_dense_2_bias_v_read_readvariableop0savev2_adam_dense_3_kernel_v_read_readvariableop.savev2_adam_dense_3_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *;
dtypes1
/2-	2
SaveV2�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes�
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

identity_1Identity_1:output:0*�
_input_shapes�
�: : : : @:@:@ : :	�@:@:@:::: : : : : : : : : : : @:@:@ : :	�@:@:@:::: : : @:@:@ : :	�@:@:@:::: 2(
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
:	�@: 
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

:: 

_output_shapes
::
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
:	�@: 

_output_shapes
:@:$ 

_output_shapes

:@: 

_output_shapes
::$ 

_output_shapes

::  

_output_shapes
::(!$
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
:	�@: (

_output_shapes
:@:$) 

_output_shapes

:@: *

_output_shapes
::$+ 

_output_shapes

:: ,

_output_shapes
::-

_output_shapes
: 
�
D
(__inference_flatten_layer_call_fn_985548

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9841782
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
e
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_984003

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim�

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+���������������������������2

ExpandDims�
MaxPoolMaxPoolExpandDims:output:0*A
_output_shapes/
-:+���������������������������*
ksize
*
paddingVALID*
strides
2	
MaxPool�
SqueezeSqueezeMaxPool:output:0*
T0*=
_output_shapes+
):'���������������������������*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'���������������������������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'���������������������������:e a
=
_output_shapes+
):'���������������������������
 
_user_specified_nameinputs
И
�

I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984609

inputs#
conv_1_984511: 
conv_1_984513: #
conv_2_984518: @
conv_2_984520:@#
conv_3_984525:@ 
conv_3_984527: !
dense_1_984533:	�@
dense_1_984535:@ 
dense_2_984538:@
dense_2_984540: 
dense_3_984543:
dense_3_984545:
identity��Conv_1/StatefulPartitionedCall�-Conv_1/bias/Regularizer/Square/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�Conv_2/StatefulPartitionedCall�-Conv_2/bias/Regularizer/Square/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�Conv_3/StatefulPartitionedCall�-Conv_3/bias/Regularizer/Square/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�Dense_1/StatefulPartitionedCall�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�Dense_2/StatefulPartitionedCall�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�Dense_3/StatefulPartitionedCall�dropout/StatefulPartitionedCall�!dropout_1/StatefulPartitionedCall�!dropout_2/StatefulPartitionedCall�
Conv_1/StatefulPartitionedCallStatefulPartitionedCallinputsconv_1_984511conv_1_984513*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:���������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9840742 
Conv_1/StatefulPartitionedCall�
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9840032
max_pooling1d/PartitionedCall�
dropout/StatefulPartitionedCallStatefulPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9844642!
dropout/StatefulPartitionedCall�
Conv_2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0conv_2_984518conv_2_984520*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������L@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9841162 
Conv_2/StatefulPartitionedCall�
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9840182!
max_pooling1d_1/PartitionedCall�
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_1/PartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9844312#
!dropout_1/StatefulPartitionedCall�
Conv_3/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0conv_3_984525conv_3_984527*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9841582 
Conv_3/StatefulPartitionedCall�
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9840332!
max_pooling1d_2/PartitionedCall�
!dropout_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling1d_2/PartitionedCall:output:0"^dropout_1/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9843982#
!dropout_2/StatefulPartitionedCall�
flatten/PartitionedCallPartitionedCall*dropout_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9841782
flatten/PartitionedCall�
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_984533dense_1_984535*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9842032!
Dense_1/StatefulPartitionedCall�
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_984538dense_2_984540*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9842322!
Dense_2/StatefulPartitionedCall�
Dense_3/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_984543dense_3_984545*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_3_layer_call_and_return_conditional_losses_9842482!
Dense_3/StatefulPartitionedCall�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984511*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984513*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984518*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984520*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984525*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984527*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984533*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984535*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984538*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984540*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentity(Dense_3/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp ^Dense_3/StatefulPartitionedCall ^dropout/StatefulPartitionedCall"^dropout_1/StatefulPartitionedCall"^dropout_2/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2@
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
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2B
Dense_3/StatefulPartitionedCallDense_3/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall2F
!dropout_2/StatefulPartitionedCall!dropout_2/StatefulPartitionedCall:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
J
.__inference_max_pooling1d_layer_call_fn_984009

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *=
_output_shapes+
):'���������������������������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9840032
PartitionedCall�
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'���������������������������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'���������������������������:e a
=
_output_shapes+
):'���������������������������
 
_user_specified_nameinputs
�
�
(__inference_Dense_3_layer_call_fn_985651

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_3_layer_call_and_return_conditional_losses_9842482
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
C__inference_Dense_1_layer_call_and_return_conditional_losses_985598

inputs1
matmul_readvariableop_resource:	�@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:���������@2
Tanh�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
F
*__inference_dropout_1_layer_call_fn_985445

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9841282
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:���������&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������&@:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�
�
(__inference_Dense_2_layer_call_fn_985619

inputs
unknown:@
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9842322
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
__inference_loss_fn_2_985694N
8conv_2_kernel_regularizer_square_readvariableop_resource: @
identity��/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp8conv_2_kernel_regularizer_square_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
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
�
�
__inference_loss_fn_5_985727D
6conv_3_bias_regularizer_square_readvariableop_resource: 
identity��-Conv_3/bias/Regularizer/Square/ReadVariableOp�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp6conv_3_bias_regularizer_square_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
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
�
d
E__inference_dropout_1_layer_call_and_return_conditional_losses_985467

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:���������&@2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:���������&@*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:���������&@2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:���������&@2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:���������&@2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:���������&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������&@:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�

�
.__inference_FCN_regressor_layer_call_fn_984991

inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	�@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_9843152
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_985315

inputsH
2conv_1_conv1d_expanddims_1_readvariableop_resource: 4
&conv_1_biasadd_readvariableop_resource: H
2conv_2_conv1d_expanddims_1_readvariableop_resource: @4
&conv_2_biasadd_readvariableop_resource:@H
2conv_3_conv1d_expanddims_1_readvariableop_resource:@ 4
&conv_3_biasadd_readvariableop_resource: 9
&dense_1_matmul_readvariableop_resource:	�@5
'dense_1_biasadd_readvariableop_resource:@8
&dense_2_matmul_readvariableop_resource:@5
'dense_2_biasadd_readvariableop_resource:8
&dense_3_matmul_readvariableop_resource:5
'dense_3_biasadd_readvariableop_resource:
identity��Conv_1/BiasAdd/ReadVariableOp�-Conv_1/bias/Regularizer/Square/ReadVariableOp�)Conv_1/conv1d/ExpandDims_1/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�Conv_2/BiasAdd/ReadVariableOp�-Conv_2/bias/Regularizer/Square/ReadVariableOp�)Conv_2/conv1d/ExpandDims_1/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�Conv_3/BiasAdd/ReadVariableOp�-Conv_3/bias/Regularizer/Square/ReadVariableOp�)Conv_3/conv1d/ExpandDims_1/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�Dense_1/BiasAdd/ReadVariableOp�Dense_1/MatMul/ReadVariableOp�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�Dense_2/BiasAdd/ReadVariableOp�Dense_2/MatMul/ReadVariableOp�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�Dense_3/BiasAdd/ReadVariableOp�Dense_3/MatMul/ReadVariableOp�
Conv_1/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
Conv_1/conv1d/ExpandDims/dim�
Conv_1/conv1d/ExpandDims
ExpandDimsinputs%Conv_1/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:����������2
Conv_1/conv1d/ExpandDims�
)Conv_1/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02+
)Conv_1/conv1d/ExpandDims_1/ReadVariableOp�
Conv_1/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_1/conv1d/ExpandDims_1/dim�
Conv_1/conv1d/ExpandDims_1
ExpandDims1Conv_1/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_1/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
Conv_1/conv1d/ExpandDims_1�
Conv_1/conv1dConv2D!Conv_1/conv1d/ExpandDims:output:0#Conv_1/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:���������� *
paddingVALID*
strides
2
Conv_1/conv1d�
Conv_1/conv1d/SqueezeSqueezeConv_1/conv1d:output:0*
T0*,
_output_shapes
:���������� *
squeeze_dims

���������2
Conv_1/conv1d/Squeeze�
Conv_1/BiasAdd/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_1/BiasAdd/ReadVariableOp�
Conv_1/BiasAddBiasAddConv_1/conv1d/Squeeze:output:0%Conv_1/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:���������� 2
Conv_1/BiasAddr
Conv_1/TanhTanhConv_1/BiasAdd:output:0*
T0*,
_output_shapes
:���������� 2
Conv_1/Tanh~
max_pooling1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
max_pooling1d/ExpandDims/dim�
max_pooling1d/ExpandDims
ExpandDimsConv_1/Tanh:y:0%max_pooling1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:���������� 2
max_pooling1d/ExpandDims�
max_pooling1d/MaxPoolMaxPool!max_pooling1d/ExpandDims:output:0*/
_output_shapes
:���������X *
ksize
*
paddingVALID*
strides
2
max_pooling1d/MaxPool�
max_pooling1d/SqueezeSqueezemax_pooling1d/MaxPool:output:0*
T0*+
_output_shapes
:���������X *
squeeze_dims
2
max_pooling1d/Squeezes
dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/dropout/Const�
dropout/dropout/MulMulmax_pooling1d/Squeeze:output:0dropout/dropout/Const:output:0*
T0*+
_output_shapes
:���������X 2
dropout/dropout/Mul|
dropout/dropout/ShapeShapemax_pooling1d/Squeeze:output:0*
T0*
_output_shapes
:2
dropout/dropout/Shape�
,dropout/dropout/random_uniform/RandomUniformRandomUniformdropout/dropout/Shape:output:0*
T0*+
_output_shapes
:���������X *
dtype02.
,dropout/dropout/random_uniform/RandomUniform�
dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2 
dropout/dropout/GreaterEqual/y�
dropout/dropout/GreaterEqualGreaterEqual5dropout/dropout/random_uniform/RandomUniform:output:0'dropout/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:���������X 2
dropout/dropout/GreaterEqual�
dropout/dropout/CastCast dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:���������X 2
dropout/dropout/Cast�
dropout/dropout/Mul_1Muldropout/dropout/Mul:z:0dropout/dropout/Cast:y:0*
T0*+
_output_shapes
:���������X 2
dropout/dropout/Mul_1�
Conv_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
Conv_2/conv1d/ExpandDims/dim�
Conv_2/conv1d/ExpandDims
ExpandDimsdropout/dropout/Mul_1:z:0%Conv_2/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������X 2
Conv_2/conv1d/ExpandDims�
)Conv_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02+
)Conv_2/conv1d/ExpandDims_1/ReadVariableOp�
Conv_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_2/conv1d/ExpandDims_1/dim�
Conv_2/conv1d/ExpandDims_1
ExpandDims1Conv_2/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_2/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
Conv_2/conv1d/ExpandDims_1�
Conv_2/conv1dConv2D!Conv_2/conv1d/ExpandDims:output:0#Conv_2/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������L@*
paddingVALID*
strides
2
Conv_2/conv1d�
Conv_2/conv1d/SqueezeSqueezeConv_2/conv1d:output:0*
T0*+
_output_shapes
:���������L@*
squeeze_dims

���������2
Conv_2/conv1d/Squeeze�
Conv_2/BiasAdd/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
Conv_2/BiasAdd/ReadVariableOp�
Conv_2/BiasAddBiasAddConv_2/conv1d/Squeeze:output:0%Conv_2/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������L@2
Conv_2/BiasAddq
Conv_2/TanhTanhConv_2/BiasAdd:output:0*
T0*+
_output_shapes
:���������L@2
Conv_2/Tanh�
max_pooling1d_1/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_1/ExpandDims/dim�
max_pooling1d_1/ExpandDims
ExpandDimsConv_2/Tanh:y:0'max_pooling1d_1/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������L@2
max_pooling1d_1/ExpandDims�
max_pooling1d_1/MaxPoolMaxPool#max_pooling1d_1/ExpandDims:output:0*/
_output_shapes
:���������&@*
ksize
*
paddingVALID*
strides
2
max_pooling1d_1/MaxPool�
max_pooling1d_1/SqueezeSqueeze max_pooling1d_1/MaxPool:output:0*
T0*+
_output_shapes
:���������&@*
squeeze_dims
2
max_pooling1d_1/Squeezew
dropout_1/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout_1/dropout/Const�
dropout_1/dropout/MulMul max_pooling1d_1/Squeeze:output:0 dropout_1/dropout/Const:output:0*
T0*+
_output_shapes
:���������&@2
dropout_1/dropout/Mul�
dropout_1/dropout/ShapeShape max_pooling1d_1/Squeeze:output:0*
T0*
_output_shapes
:2
dropout_1/dropout/Shape�
.dropout_1/dropout/random_uniform/RandomUniformRandomUniform dropout_1/dropout/Shape:output:0*
T0*+
_output_shapes
:���������&@*
dtype020
.dropout_1/dropout/random_uniform/RandomUniform�
 dropout_1/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2"
 dropout_1/dropout/GreaterEqual/y�
dropout_1/dropout/GreaterEqualGreaterEqual7dropout_1/dropout/random_uniform/RandomUniform:output:0)dropout_1/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:���������&@2 
dropout_1/dropout/GreaterEqual�
dropout_1/dropout/CastCast"dropout_1/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:���������&@2
dropout_1/dropout/Cast�
dropout_1/dropout/Mul_1Muldropout_1/dropout/Mul:z:0dropout_1/dropout/Cast:y:0*
T0*+
_output_shapes
:���������&@2
dropout_1/dropout/Mul_1�
Conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
Conv_3/conv1d/ExpandDims/dim�
Conv_3/conv1d/ExpandDims
ExpandDimsdropout_1/dropout/Mul_1:z:0%Conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������&@2
Conv_3/conv1d/ExpandDims�
)Conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02+
)Conv_3/conv1d/ExpandDims_1/ReadVariableOp�
Conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_3/conv1d/ExpandDims_1/dim�
Conv_3/conv1d/ExpandDims_1
ExpandDims1Conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
Conv_3/conv1d/ExpandDims_1�
Conv_3/conv1dConv2D!Conv_3/conv1d/ExpandDims:output:0#Conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������  *
paddingVALID*
strides
2
Conv_3/conv1d�
Conv_3/conv1d/SqueezeSqueezeConv_3/conv1d:output:0*
T0*+
_output_shapes
:���������  *
squeeze_dims

���������2
Conv_3/conv1d/Squeeze�
Conv_3/BiasAdd/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_3/BiasAdd/ReadVariableOp�
Conv_3/BiasAddBiasAddConv_3/conv1d/Squeeze:output:0%Conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������  2
Conv_3/BiasAddq
Conv_3/TanhTanhConv_3/BiasAdd:output:0*
T0*+
_output_shapes
:���������  2
Conv_3/Tanh�
max_pooling1d_2/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_2/ExpandDims/dim�
max_pooling1d_2/ExpandDims
ExpandDimsConv_3/Tanh:y:0'max_pooling1d_2/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������  2
max_pooling1d_2/ExpandDims�
max_pooling1d_2/MaxPoolMaxPool#max_pooling1d_2/ExpandDims:output:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2
max_pooling1d_2/MaxPool�
max_pooling1d_2/SqueezeSqueeze max_pooling1d_2/MaxPool:output:0*
T0*+
_output_shapes
:��������� *
squeeze_dims
2
max_pooling1d_2/Squeezew
dropout_2/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout_2/dropout/Const�
dropout_2/dropout/MulMul max_pooling1d_2/Squeeze:output:0 dropout_2/dropout/Const:output:0*
T0*+
_output_shapes
:��������� 2
dropout_2/dropout/Mul�
dropout_2/dropout/ShapeShape max_pooling1d_2/Squeeze:output:0*
T0*
_output_shapes
:2
dropout_2/dropout/Shape�
.dropout_2/dropout/random_uniform/RandomUniformRandomUniform dropout_2/dropout/Shape:output:0*
T0*+
_output_shapes
:��������� *
dtype020
.dropout_2/dropout/random_uniform/RandomUniform�
 dropout_2/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2"
 dropout_2/dropout/GreaterEqual/y�
dropout_2/dropout/GreaterEqualGreaterEqual7dropout_2/dropout/random_uniform/RandomUniform:output:0)dropout_2/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:��������� 2 
dropout_2/dropout/GreaterEqual�
dropout_2/dropout/CastCast"dropout_2/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:��������� 2
dropout_2/dropout/Cast�
dropout_2/dropout/Mul_1Muldropout_2/dropout/Mul:z:0dropout_2/dropout/Cast:y:0*
T0*+
_output_shapes
:��������� 2
dropout_2/dropout/Mul_1o
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����   2
flatten/Const�
flatten/ReshapeReshapedropout_2/dropout/Mul_1:z:0flatten/Const:output:0*
T0*(
_output_shapes
:����������2
flatten/Reshape�
Dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	�@*
dtype02
Dense_1/MatMul/ReadVariableOp�
Dense_1/MatMulMatMulflatten/Reshape:output:0%Dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
Dense_1/MatMul�
Dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02 
Dense_1/BiasAdd/ReadVariableOp�
Dense_1/BiasAddBiasAddDense_1/MatMul:product:0&Dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
Dense_1/BiasAddp
Dense_1/TanhTanhDense_1/BiasAdd:output:0*
T0*'
_output_shapes
:���������@2
Dense_1/Tanh�
Dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype02
Dense_2/MatMul/ReadVariableOp�
Dense_2/MatMulMatMulDense_1/Tanh:y:0%Dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_2/MatMul�
Dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
Dense_2/BiasAdd/ReadVariableOp�
Dense_2/BiasAddBiasAddDense_2/MatMul:product:0&Dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_2/BiasAddp
Dense_2/TanhTanhDense_2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
Dense_2/Tanh�
Dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
Dense_3/MatMul/ReadVariableOp�
Dense_3/MatMulMatMulDense_2/Tanh:y:0%Dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_3/MatMul�
Dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
Dense_3/BiasAdd/ReadVariableOp�
Dense_3/BiasAddBiasAddDense_3/MatMul:product:0&Dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_3/BiasAdd�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentityDense_3/BiasAdd:output:0^Conv_1/BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp*^Conv_1/conv1d/ExpandDims_1/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp*^Conv_2/conv1d/ExpandDims_1/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp*^Conv_3/conv1d/ExpandDims_1/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp^Dense_1/BiasAdd/ReadVariableOp^Dense_1/MatMul/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^Dense_2/BiasAdd/ReadVariableOp^Dense_2/MatMul/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp^Dense_3/BiasAdd/ReadVariableOp^Dense_3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2>
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
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2@
Dense_3/BiasAdd/ReadVariableOpDense_3/BiasAdd/ReadVariableOp2>
Dense_3/MatMul/ReadVariableOpDense_3/MatMul/ReadVariableOp:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
F
*__inference_dropout_2_layer_call_fn_985521

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9841702
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:��������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
�
__inference_loss_fn_6_985738L
9dense_1_kernel_regularizer_square_readvariableop_resource:	�@
identity��0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_1_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
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
�
�
'__inference_Conv_3_layer_call_fn_985488

inputs
unknown:@ 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9841582
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:���������  2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������&@: : 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�
L
0__inference_max_pooling1d_1_layer_call_fn_984024

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *=
_output_shapes+
):'���������������������������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9840182
PartitionedCall�
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'���������������������������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*<
_input_shapes+
):'���������������������������:e a
=
_output_shapes+
):'���������������������������
 
_user_specified_nameinputs
�
�
C__inference_Dense_1_layer_call_and_return_conditional_losses_984203

inputs1
matmul_readvariableop_resource:	�@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:���������@2
Tanh�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp2d
0Dense_1/kernel/Regularizer/Square/ReadVariableOp0Dense_1/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
d
E__inference_dropout_2_layer_call_and_return_conditional_losses_985543

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:��������� 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:��������� *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:��������� 2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:��������� 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:��������� 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:��������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
c
E__inference_dropout_2_layer_call_and_return_conditional_losses_984170

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:��������� 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:��������� 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
c
E__inference_dropout_1_layer_call_and_return_conditional_losses_985455

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:���������&@2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:���������&@2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������&@:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�

�
$__inference_signature_wrapper_984962
convolutional_inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	�@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__wrapped_model_9839942
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
,
_output_shapes
:����������
.
_user_specified_nameConvolutional_inputs
�
c
E__inference_dropout_2_layer_call_and_return_conditional_losses_985531

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:��������� 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:��������� 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� :S O
+
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
�
__inference_loss_fn_0_985672N
8conv_1_kernel_regularizer_square_readvariableop_resource: 
identity��/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp8conv_1_kernel_regularizer_square_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
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
�
�
(__inference_Dense_1_layer_call_fn_985575

inputs
unknown:	�@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9842032
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�	
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984766
convolutional_inputs#
conv_1_984668: 
conv_1_984670: #
conv_2_984675: @
conv_2_984677:@#
conv_3_984682:@ 
conv_3_984684: !
dense_1_984690:	�@
dense_1_984692:@ 
dense_2_984695:@
dense_2_984697: 
dense_3_984700:
dense_3_984702:
identity��Conv_1/StatefulPartitionedCall�-Conv_1/bias/Regularizer/Square/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�Conv_2/StatefulPartitionedCall�-Conv_2/bias/Regularizer/Square/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�Conv_3/StatefulPartitionedCall�-Conv_3/bias/Regularizer/Square/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�Dense_1/StatefulPartitionedCall�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�Dense_2/StatefulPartitionedCall�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�Dense_3/StatefulPartitionedCall�
Conv_1/StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsconv_1_984668conv_1_984670*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:���������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_1_layer_call_and_return_conditional_losses_9840742 
Conv_1/StatefulPartitionedCall�
max_pooling1d/PartitionedCallPartitionedCall'Conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_9840032
max_pooling1d/PartitionedCall�
dropout/PartitionedCallPartitionedCall&max_pooling1d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������X * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_dropout_layer_call_and_return_conditional_losses_9840862
dropout/PartitionedCall�
Conv_2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0conv_2_984675conv_2_984677*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������L@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9841162 
Conv_2/StatefulPartitionedCall�
max_pooling1d_1/PartitionedCallPartitionedCall'Conv_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_9840182!
max_pooling1d_1/PartitionedCall�
dropout_1/PartitionedCallPartitionedCall(max_pooling1d_1/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������&@* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_1_layer_call_and_return_conditional_losses_9841282
dropout_1/PartitionedCall�
Conv_3/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0conv_3_984682conv_3_984684*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������  *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_3_layer_call_and_return_conditional_losses_9841582 
Conv_3/StatefulPartitionedCall�
max_pooling1d_2/PartitionedCallPartitionedCall'Conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_9840332!
max_pooling1d_2/PartitionedCall�
dropout_2/PartitionedCallPartitionedCall(max_pooling1d_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dropout_2_layer_call_and_return_conditional_losses_9841702
dropout_2/PartitionedCall�
flatten/PartitionedCallPartitionedCall"dropout_2/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_flatten_layer_call_and_return_conditional_losses_9841782
flatten/PartitionedCall�
Dense_1/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_1_984690dense_1_984692*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_1_layer_call_and_return_conditional_losses_9842032!
Dense_1/StatefulPartitionedCall�
Dense_2/StatefulPartitionedCallStatefulPartitionedCall(Dense_1/StatefulPartitionedCall:output:0dense_2_984695dense_2_984697*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_2_layer_call_and_return_conditional_losses_9842322!
Dense_2/StatefulPartitionedCall�
Dense_3/StatefulPartitionedCallStatefulPartitionedCall(Dense_2/StatefulPartitionedCall:output:0dense_3_984700dense_3_984702*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_Dense_3_layer_call_and_return_conditional_losses_9842482!
Dense_3/StatefulPartitionedCall�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984668*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_1_984670*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984675*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_2_984677*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984682*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv_3_984684*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984690*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_984692*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984695*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_984697*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentity(Dense_3/StatefulPartitionedCall:output:0^Conv_1/StatefulPartitionedCall.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/StatefulPartitionedCall.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/StatefulPartitionedCall.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp ^Dense_1/StatefulPartitionedCall/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp ^Dense_2/StatefulPartitionedCall/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp ^Dense_3/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2@
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
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2B
Dense_3/StatefulPartitionedCallDense_3/StatefulPartitionedCall:b ^
,
_output_shapes
:����������
.
_user_specified_nameConvolutional_inputs
�
�
.__inference_FCN_regressor_layer_call_fn_984665
convolutional_inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	�@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_9846092
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
,
_output_shapes
:����������
.
_user_specified_nameConvolutional_inputs
�
b
C__inference_dropout_layer_call_and_return_conditional_losses_985391

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:���������X 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:���������X *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:���������X 2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:���������X 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:���������X 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:���������X 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������X :S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�
�
'__inference_Conv_2_layer_call_fn_985412

inputs
unknown: @
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:���������L@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_Conv_2_layer_call_and_return_conditional_losses_9841162
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:���������L@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������X : : 22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�$
�
B__inference_Conv_2_layer_call_and_return_conditional_losses_984116

inputsA
+conv1d_expanddims_1_readvariableop_resource: @-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�-Conv_2/bias/Regularizer/Square/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
conv1d/ExpandDims/dim�
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������X 2
conv1d/ExpandDims�
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
conv1d/ExpandDims_1/dim�
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
conv1d/ExpandDims_1�
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������L@*
paddingVALID*
strides
2
conv1d�
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:���������L@*
squeeze_dims

���������2
conv1d/Squeeze�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������L@2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:���������L@2
Tanh�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:���������L@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������X : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�
b
C__inference_dropout_layer_call_and_return_conditional_losses_984464

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:���������X 2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:���������X *
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:���������X 2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:���������X 2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:���������X 2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:���������X 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������X :S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
��
�
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_985157

inputsH
2conv_1_conv1d_expanddims_1_readvariableop_resource: 4
&conv_1_biasadd_readvariableop_resource: H
2conv_2_conv1d_expanddims_1_readvariableop_resource: @4
&conv_2_biasadd_readvariableop_resource:@H
2conv_3_conv1d_expanddims_1_readvariableop_resource:@ 4
&conv_3_biasadd_readvariableop_resource: 9
&dense_1_matmul_readvariableop_resource:	�@5
'dense_1_biasadd_readvariableop_resource:@8
&dense_2_matmul_readvariableop_resource:@5
'dense_2_biasadd_readvariableop_resource:8
&dense_3_matmul_readvariableop_resource:5
'dense_3_biasadd_readvariableop_resource:
identity��Conv_1/BiasAdd/ReadVariableOp�-Conv_1/bias/Regularizer/Square/ReadVariableOp�)Conv_1/conv1d/ExpandDims_1/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�Conv_2/BiasAdd/ReadVariableOp�-Conv_2/bias/Regularizer/Square/ReadVariableOp�)Conv_2/conv1d/ExpandDims_1/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�Conv_3/BiasAdd/ReadVariableOp�-Conv_3/bias/Regularizer/Square/ReadVariableOp�)Conv_3/conv1d/ExpandDims_1/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�Dense_1/BiasAdd/ReadVariableOp�Dense_1/MatMul/ReadVariableOp�.Dense_1/bias/Regularizer/Square/ReadVariableOp�0Dense_1/kernel/Regularizer/Square/ReadVariableOp�Dense_2/BiasAdd/ReadVariableOp�Dense_2/MatMul/ReadVariableOp�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�Dense_3/BiasAdd/ReadVariableOp�Dense_3/MatMul/ReadVariableOp�
Conv_1/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
Conv_1/conv1d/ExpandDims/dim�
Conv_1/conv1d/ExpandDims
ExpandDimsinputs%Conv_1/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:����������2
Conv_1/conv1d/ExpandDims�
)Conv_1/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype02+
)Conv_1/conv1d/ExpandDims_1/ReadVariableOp�
Conv_1/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_1/conv1d/ExpandDims_1/dim�
Conv_1/conv1d/ExpandDims_1
ExpandDims1Conv_1/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_1/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
Conv_1/conv1d/ExpandDims_1�
Conv_1/conv1dConv2D!Conv_1/conv1d/ExpandDims:output:0#Conv_1/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:���������� *
paddingVALID*
strides
2
Conv_1/conv1d�
Conv_1/conv1d/SqueezeSqueezeConv_1/conv1d:output:0*
T0*,
_output_shapes
:���������� *
squeeze_dims

���������2
Conv_1/conv1d/Squeeze�
Conv_1/BiasAdd/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_1/BiasAdd/ReadVariableOp�
Conv_1/BiasAddBiasAddConv_1/conv1d/Squeeze:output:0%Conv_1/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:���������� 2
Conv_1/BiasAddr
Conv_1/TanhTanhConv_1/BiasAdd:output:0*
T0*,
_output_shapes
:���������� 2
Conv_1/Tanh~
max_pooling1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
max_pooling1d/ExpandDims/dim�
max_pooling1d/ExpandDims
ExpandDimsConv_1/Tanh:y:0%max_pooling1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:���������� 2
max_pooling1d/ExpandDims�
max_pooling1d/MaxPoolMaxPool!max_pooling1d/ExpandDims:output:0*/
_output_shapes
:���������X *
ksize
*
paddingVALID*
strides
2
max_pooling1d/MaxPool�
max_pooling1d/SqueezeSqueezemax_pooling1d/MaxPool:output:0*
T0*+
_output_shapes
:���������X *
squeeze_dims
2
max_pooling1d/Squeeze�
dropout/IdentityIdentitymax_pooling1d/Squeeze:output:0*
T0*+
_output_shapes
:���������X 2
dropout/Identity�
Conv_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
Conv_2/conv1d/ExpandDims/dim�
Conv_2/conv1d/ExpandDims
ExpandDimsdropout/Identity:output:0%Conv_2/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������X 2
Conv_2/conv1d/ExpandDims�
)Conv_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype02+
)Conv_2/conv1d/ExpandDims_1/ReadVariableOp�
Conv_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_2/conv1d/ExpandDims_1/dim�
Conv_2/conv1d/ExpandDims_1
ExpandDims1Conv_2/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_2/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
Conv_2/conv1d/ExpandDims_1�
Conv_2/conv1dConv2D!Conv_2/conv1d/ExpandDims:output:0#Conv_2/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������L@*
paddingVALID*
strides
2
Conv_2/conv1d�
Conv_2/conv1d/SqueezeSqueezeConv_2/conv1d:output:0*
T0*+
_output_shapes
:���������L@*
squeeze_dims

���������2
Conv_2/conv1d/Squeeze�
Conv_2/BiasAdd/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
Conv_2/BiasAdd/ReadVariableOp�
Conv_2/BiasAddBiasAddConv_2/conv1d/Squeeze:output:0%Conv_2/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������L@2
Conv_2/BiasAddq
Conv_2/TanhTanhConv_2/BiasAdd:output:0*
T0*+
_output_shapes
:���������L@2
Conv_2/Tanh�
max_pooling1d_1/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_1/ExpandDims/dim�
max_pooling1d_1/ExpandDims
ExpandDimsConv_2/Tanh:y:0'max_pooling1d_1/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������L@2
max_pooling1d_1/ExpandDims�
max_pooling1d_1/MaxPoolMaxPool#max_pooling1d_1/ExpandDims:output:0*/
_output_shapes
:���������&@*
ksize
*
paddingVALID*
strides
2
max_pooling1d_1/MaxPool�
max_pooling1d_1/SqueezeSqueeze max_pooling1d_1/MaxPool:output:0*
T0*+
_output_shapes
:���������&@*
squeeze_dims
2
max_pooling1d_1/Squeeze�
dropout_1/IdentityIdentity max_pooling1d_1/Squeeze:output:0*
T0*+
_output_shapes
:���������&@2
dropout_1/Identity�
Conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
Conv_3/conv1d/ExpandDims/dim�
Conv_3/conv1d/ExpandDims
ExpandDimsdropout_1/Identity:output:0%Conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������&@2
Conv_3/conv1d/ExpandDims�
)Conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype02+
)Conv_3/conv1d/ExpandDims_1/ReadVariableOp�
Conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
Conv_3/conv1d/ExpandDims_1/dim�
Conv_3/conv1d/ExpandDims_1
ExpandDims1Conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0'Conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
Conv_3/conv1d/ExpandDims_1�
Conv_3/conv1dConv2D!Conv_3/conv1d/ExpandDims:output:0#Conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������  *
paddingVALID*
strides
2
Conv_3/conv1d�
Conv_3/conv1d/SqueezeSqueezeConv_3/conv1d:output:0*
T0*+
_output_shapes
:���������  *
squeeze_dims

���������2
Conv_3/conv1d/Squeeze�
Conv_3/BiasAdd/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02
Conv_3/BiasAdd/ReadVariableOp�
Conv_3/BiasAddBiasAddConv_3/conv1d/Squeeze:output:0%Conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������  2
Conv_3/BiasAddq
Conv_3/TanhTanhConv_3/BiasAdd:output:0*
T0*+
_output_shapes
:���������  2
Conv_3/Tanh�
max_pooling1d_2/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2 
max_pooling1d_2/ExpandDims/dim�
max_pooling1d_2/ExpandDims
ExpandDimsConv_3/Tanh:y:0'max_pooling1d_2/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������  2
max_pooling1d_2/ExpandDims�
max_pooling1d_2/MaxPoolMaxPool#max_pooling1d_2/ExpandDims:output:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2
max_pooling1d_2/MaxPool�
max_pooling1d_2/SqueezeSqueeze max_pooling1d_2/MaxPool:output:0*
T0*+
_output_shapes
:��������� *
squeeze_dims
2
max_pooling1d_2/Squeeze�
dropout_2/IdentityIdentity max_pooling1d_2/Squeeze:output:0*
T0*+
_output_shapes
:��������� 2
dropout_2/Identityo
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����   2
flatten/Const�
flatten/ReshapeReshapedropout_2/Identity:output:0flatten/Const:output:0*
T0*(
_output_shapes
:����������2
flatten/Reshape�
Dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	�@*
dtype02
Dense_1/MatMul/ReadVariableOp�
Dense_1/MatMulMatMulflatten/Reshape:output:0%Dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
Dense_1/MatMul�
Dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02 
Dense_1/BiasAdd/ReadVariableOp�
Dense_1/BiasAddBiasAddDense_1/MatMul:product:0&Dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@2
Dense_1/BiasAddp
Dense_1/TanhTanhDense_1/BiasAdd:output:0*
T0*'
_output_shapes
:���������@2
Dense_1/Tanh�
Dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype02
Dense_2/MatMul/ReadVariableOp�
Dense_2/MatMulMatMulDense_1/Tanh:y:0%Dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_2/MatMul�
Dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
Dense_2/BiasAdd/ReadVariableOp�
Dense_2/BiasAddBiasAddDense_2/MatMul:product:0&Dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_2/BiasAddp
Dense_2/TanhTanhDense_2/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
Dense_2/Tanh�
Dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
Dense_3/MatMul/ReadVariableOp�
Dense_3/MatMulMatMulDense_2/Tanh:y:0%Dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_3/MatMul�
Dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
Dense_3/BiasAdd/ReadVariableOp�
Dense_3/BiasAddBiasAddDense_3/MatMul:product:0&Dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
Dense_3/BiasAdd�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_1_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_2_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
0Dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	�@*
dtype022
0Dense_1/kernel/Regularizer/Square/ReadVariableOp�
!Dense_1/kernel/Regularizer/SquareSquare8Dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�@2#
!Dense_1/kernel/Regularizer/Square�
 Dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_1/kernel/Regularizer/Const�
Dense_1/kernel/Regularizer/SumSum%Dense_1/kernel/Regularizer/Square:y:0)Dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/Sum�
 Dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_1/kernel/Regularizer/mul/x�
Dense_1/kernel/Regularizer/mulMul)Dense_1/kernel/Regularizer/mul/x:output:0'Dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_1/kernel/Regularizer/mul�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentityDense_3/BiasAdd:output:0^Conv_1/BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp*^Conv_1/conv1d/ExpandDims_1/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp^Conv_2/BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp*^Conv_2/conv1d/ExpandDims_1/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp^Conv_3/BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp*^Conv_3/conv1d/ExpandDims_1/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp^Dense_1/BiasAdd/ReadVariableOp^Dense_1/MatMul/ReadVariableOp/^Dense_1/bias/Regularizer/Square/ReadVariableOp1^Dense_1/kernel/Regularizer/Square/ReadVariableOp^Dense_2/BiasAdd/ReadVariableOp^Dense_2/MatMul/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp^Dense_3/BiasAdd/ReadVariableOp^Dense_3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 2>
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
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2@
Dense_3/BiasAdd/ReadVariableOpDense_3/BiasAdd/ReadVariableOp2>
Dense_3/MatMul/ReadVariableOpDense_3/MatMul/ReadVariableOp:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
a
C__inference_dropout_layer_call_and_return_conditional_losses_985379

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:���������X 2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:���������X 2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������X :S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�$
�
B__inference_Conv_1_layer_call_and_return_conditional_losses_985364

inputsA
+conv1d_expanddims_1_readvariableop_resource: -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�-Conv_1/bias/Regularizer/Square/ReadVariableOp�/Conv_1/kernel/Regularizer/Square/ReadVariableOp�"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
conv1d/ExpandDims/dim�
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:����������2
conv1d/ExpandDims�
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
conv1d/ExpandDims_1/dim�
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: 2
conv1d/ExpandDims_1�
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:���������� *
paddingVALID*
strides
2
conv1d�
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:���������� *
squeeze_dims

���������2
conv1d/Squeeze�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:���������� 2	
BiasAdd]
TanhTanhBiasAdd:output:0*
T0*,
_output_shapes
:���������� 2
Tanh�
/Conv_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: *
dtype021
/Conv_1/kernel/Regularizer/Square/ReadVariableOp�
 Conv_1/kernel/Regularizer/SquareSquare7Conv_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: 2"
 Conv_1/kernel/Regularizer/Square�
Conv_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_1/kernel/Regularizer/Const�
Conv_1/kernel/Regularizer/SumSum$Conv_1/kernel/Regularizer/Square:y:0(Conv_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/Sum�
Conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_1/kernel/Regularizer/mul/x�
Conv_1/kernel/Regularizer/mulMul(Conv_1/kernel/Regularizer/mul/x:output:0&Conv_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/kernel/Regularizer/mul�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_1/bias/Regularizer/Square/ReadVariableOp0^Conv_1/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:���������� 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_1/bias/Regularizer/Square/ReadVariableOp-Conv_1/bias/Regularizer/Square/ReadVariableOp2b
/Conv_1/kernel/Regularizer/Square/ReadVariableOp/Conv_1/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_3_985705D
6conv_2_bias_regularizer_square_readvariableop_resource:@
identity��-Conv_2/bias/Regularizer/Square/ReadVariableOp�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp6conv_2_bias_regularizer_square_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
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
�$
�
B__inference_Conv_2_layer_call_and_return_conditional_losses_985440

inputsA
+conv1d_expanddims_1_readvariableop_resource: @-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�-Conv_2/bias/Regularizer/Square/ReadVariableOp�/Conv_2/kernel/Regularizer/Square/ReadVariableOp�"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
conv1d/ExpandDims/dim�
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������X 2
conv1d/ExpandDims�
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
conv1d/ExpandDims_1/dim�
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
: @2
conv1d/ExpandDims_1�
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������L@*
paddingVALID*
strides
2
conv1d�
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:���������L@*
squeeze_dims

���������2
conv1d/Squeeze�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������L@2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:���������L@2
Tanh�
/Conv_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
: @*
dtype021
/Conv_2/kernel/Regularizer/Square/ReadVariableOp�
 Conv_2/kernel/Regularizer/SquareSquare7Conv_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
: @2"
 Conv_2/kernel/Regularizer/Square�
Conv_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_2/kernel/Regularizer/Const�
Conv_2/kernel/Regularizer/SumSum$Conv_2/kernel/Regularizer/Square:y:0(Conv_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/Sum�
Conv_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_2/kernel/Regularizer/mul/x�
Conv_2/kernel/Regularizer/mulMul(Conv_2/kernel/Regularizer/mul/x:output:0&Conv_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/kernel/Regularizer/mul�
-Conv_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02/
-Conv_2/bias/Regularizer/Square/ReadVariableOp�
Conv_2/bias/Regularizer/SquareSquare5Conv_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2 
Conv_2/bias/Regularizer/Square�
Conv_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_2/bias/Regularizer/Const�
Conv_2/bias/Regularizer/SumSum"Conv_2/bias/Regularizer/Square:y:0&Conv_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/Sum�
Conv_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_2/bias/Regularizer/mul/x�
Conv_2/bias/Regularizer/mulMul&Conv_2/bias/Regularizer/mul/x:output:0$Conv_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_2/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_2/bias/Regularizer/Square/ReadVariableOp0^Conv_2/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:���������L@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������X : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_2/bias/Regularizer/Square/ReadVariableOp-Conv_2/bias/Regularizer/Square/ReadVariableOp2b
/Conv_2/kernel/Regularizer/Square/ReadVariableOp/Conv_2/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:���������X 
 
_user_specified_nameinputs
�
�
__inference_loss_fn_9_985771E
7dense_2_bias_regularizer_square_readvariableop_resource:
identity��.Dense_2/bias/Regularizer/Square/ReadVariableOp�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp7dense_2_bias_regularizer_square_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
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
�
�
.__inference_FCN_regressor_layer_call_fn_984342
convolutional_inputs
unknown: 
	unknown_0: 
	unknown_1: @
	unknown_2:@
	unknown_3:@ 
	unknown_4: 
	unknown_5:	�@
	unknown_6:@
	unknown_7:@
	unknown_8:
	unknown_9:

unknown_10:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallconvolutional_inputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_9843152
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
,
_output_shapes
:����������
.
_user_specified_nameConvolutional_inputs
�
�
C__inference_Dense_2_layer_call_and_return_conditional_losses_984232

inputs0
matmul_readvariableop_resource:@-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�.Dense_2/bias/Regularizer/Square/ReadVariableOp�0Dense_2/kernel/Regularizer/Square/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Tanh�
0Dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype022
0Dense_2/kernel/Regularizer/Square/ReadVariableOp�
!Dense_2/kernel/Regularizer/SquareSquare8Dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:@2#
!Dense_2/kernel/Regularizer/Square�
 Dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 Dense_2/kernel/Regularizer/Const�
Dense_2/kernel/Regularizer/SumSum%Dense_2/kernel/Regularizer/Square:y:0)Dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/Sum�
 Dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 Dense_2/kernel/Regularizer/mul/x�
Dense_2/kernel/Regularizer/mulMul)Dense_2/kernel/Regularizer/mul/x:output:0'Dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
Dense_2/kernel/Regularizer/mul�
.Dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.Dense_2/bias/Regularizer/Square/ReadVariableOp�
Dense_2/bias/Regularizer/SquareSquare6Dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:2!
Dense_2/bias/Regularizer/Square�
Dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_2/bias/Regularizer/Const�
Dense_2/bias/Regularizer/SumSum#Dense_2/bias/Regularizer/Square:y:0'Dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/Sum�
Dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_2/bias/Regularizer/mul/x�
Dense_2/bias/Regularizer/mulMul'Dense_2/bias/Regularizer/mul/x:output:0%Dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_2/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp/^Dense_2/bias/Regularizer/Square/ReadVariableOp1^Dense_2/kernel/Regularizer/Square/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2`
.Dense_2/bias/Regularizer/Square/ReadVariableOp.Dense_2/bias/Regularizer/Square/ReadVariableOp2d
0Dense_2/kernel/Regularizer/Square/ReadVariableOp0Dense_2/kernel/Regularizer/Square/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�$
�
B__inference_Conv_3_layer_call_and_return_conditional_losses_984158

inputsA
+conv1d_expanddims_1_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�-Conv_3/bias/Regularizer/Square/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
conv1d/ExpandDims/dim�
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������&@2
conv1d/ExpandDims�
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
conv1d/ExpandDims_1/dim�
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
conv1d/ExpandDims_1�
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������  *
paddingVALID*
strides
2
conv1d�
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:���������  *
squeeze_dims

���������2
conv1d/Squeeze�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������  2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:���������  2
Tanh�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:���������  2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������&@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�
d
E__inference_dropout_1_layer_call_and_return_conditional_losses_984431

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:���������&@2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:���������&@*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:���������&@2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:���������&@2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:���������&@2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:���������&@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������&@:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�
�
"__inference__traced_restore_986068
file_prefix4
assignvariableop_conv_1_kernel: ,
assignvariableop_1_conv_1_bias: 6
 assignvariableop_2_conv_2_kernel: @,
assignvariableop_3_conv_2_bias:@6
 assignvariableop_4_conv_3_kernel:@ ,
assignvariableop_5_conv_3_bias: 4
!assignvariableop_6_dense_1_kernel:	�@-
assignvariableop_7_dense_1_bias:@3
!assignvariableop_8_dense_2_kernel:@-
assignvariableop_9_dense_2_bias:4
"assignvariableop_10_dense_3_kernel:.
 assignvariableop_11_dense_3_bias:'
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
)assignvariableop_26_adam_dense_1_kernel_m:	�@5
'assignvariableop_27_adam_dense_1_bias_m:@;
)assignvariableop_28_adam_dense_2_kernel_m:@5
'assignvariableop_29_adam_dense_2_bias_m:;
)assignvariableop_30_adam_dense_3_kernel_m:5
'assignvariableop_31_adam_dense_3_bias_m:>
(assignvariableop_32_adam_conv_1_kernel_v: 4
&assignvariableop_33_adam_conv_1_bias_v: >
(assignvariableop_34_adam_conv_2_kernel_v: @4
&assignvariableop_35_adam_conv_2_bias_v:@>
(assignvariableop_36_adam_conv_3_kernel_v:@ 4
&assignvariableop_37_adam_conv_3_bias_v: <
)assignvariableop_38_adam_dense_1_kernel_v:	�@5
'assignvariableop_39_adam_dense_1_bias_v:@;
)assignvariableop_40_adam_dense_2_kernel_v:@5
'assignvariableop_41_adam_dense_2_bias_v:;
)assignvariableop_42_adam_dense_3_kernel_v:5
'assignvariableop_43_adam_dense_3_bias_v:
identity_45��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*�
value�B�-B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:-*
dtype0*m
valuedBb-B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices�
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�:::::::::::::::::::::::::::::::::::::::::::::*;
dtypes1
/2-	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

Identity�
AssignVariableOpAssignVariableOpassignvariableop_conv_1_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1�
AssignVariableOp_1AssignVariableOpassignvariableop_1_conv_1_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2�
AssignVariableOp_2AssignVariableOp assignvariableop_2_conv_2_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3�
AssignVariableOp_3AssignVariableOpassignvariableop_3_conv_2_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4�
AssignVariableOp_4AssignVariableOp assignvariableop_4_conv_3_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5�
AssignVariableOp_5AssignVariableOpassignvariableop_5_conv_3_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6�
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_1_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7�
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_1_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8�
AssignVariableOp_8AssignVariableOp!assignvariableop_8_dense_2_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9�
AssignVariableOp_9AssignVariableOpassignvariableop_9_dense_2_biasIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10�
AssignVariableOp_10AssignVariableOp"assignvariableop_10_dense_3_kernelIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11�
AssignVariableOp_11AssignVariableOp assignvariableop_11_dense_3_biasIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_12�
AssignVariableOp_12AssignVariableOpassignvariableop_12_adam_iterIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13�
AssignVariableOp_13AssignVariableOpassignvariableop_13_adam_beta_1Identity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14�
AssignVariableOp_14AssignVariableOpassignvariableop_14_adam_beta_2Identity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15�
AssignVariableOp_15AssignVariableOpassignvariableop_15_adam_decayIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16�
AssignVariableOp_16AssignVariableOpassignvariableop_16_totalIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17�
AssignVariableOp_17AssignVariableOpassignvariableop_17_countIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18�
AssignVariableOp_18AssignVariableOpassignvariableop_18_total_1Identity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19�
AssignVariableOp_19AssignVariableOpassignvariableop_19_count_1Identity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20�
AssignVariableOp_20AssignVariableOp(assignvariableop_20_adam_conv_1_kernel_mIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21�
AssignVariableOp_21AssignVariableOp&assignvariableop_21_adam_conv_1_bias_mIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22�
AssignVariableOp_22AssignVariableOp(assignvariableop_22_adam_conv_2_kernel_mIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23�
AssignVariableOp_23AssignVariableOp&assignvariableop_23_adam_conv_2_bias_mIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24�
AssignVariableOp_24AssignVariableOp(assignvariableop_24_adam_conv_3_kernel_mIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25�
AssignVariableOp_25AssignVariableOp&assignvariableop_25_adam_conv_3_bias_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26�
AssignVariableOp_26AssignVariableOp)assignvariableop_26_adam_dense_1_kernel_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27�
AssignVariableOp_27AssignVariableOp'assignvariableop_27_adam_dense_1_bias_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28�
AssignVariableOp_28AssignVariableOp)assignvariableop_28_adam_dense_2_kernel_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29�
AssignVariableOp_29AssignVariableOp'assignvariableop_29_adam_dense_2_bias_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30�
AssignVariableOp_30AssignVariableOp)assignvariableop_30_adam_dense_3_kernel_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31�
AssignVariableOp_31AssignVariableOp'assignvariableop_31_adam_dense_3_bias_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32�
AssignVariableOp_32AssignVariableOp(assignvariableop_32_adam_conv_1_kernel_vIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33�
AssignVariableOp_33AssignVariableOp&assignvariableop_33_adam_conv_1_bias_vIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34�
AssignVariableOp_34AssignVariableOp(assignvariableop_34_adam_conv_2_kernel_vIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35�
AssignVariableOp_35AssignVariableOp&assignvariableop_35_adam_conv_2_bias_vIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36�
AssignVariableOp_36AssignVariableOp(assignvariableop_36_adam_conv_3_kernel_vIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37�
AssignVariableOp_37AssignVariableOp&assignvariableop_37_adam_conv_3_bias_vIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38�
AssignVariableOp_38AssignVariableOp)assignvariableop_38_adam_dense_1_kernel_vIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39�
AssignVariableOp_39AssignVariableOp'assignvariableop_39_adam_dense_1_bias_vIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40�
AssignVariableOp_40AssignVariableOp)assignvariableop_40_adam_dense_2_kernel_vIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41�
AssignVariableOp_41AssignVariableOp'assignvariableop_41_adam_dense_2_bias_vIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42�
AssignVariableOp_42AssignVariableOp)assignvariableop_42_adam_dense_3_kernel_vIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43�
AssignVariableOp_43AssignVariableOp'assignvariableop_43_adam_dense_3_bias_vIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_439
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp�
Identity_44Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_44�
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
�
�
__inference_loss_fn_1_985683D
6conv_1_bias_regularizer_square_readvariableop_resource: 
identity��-Conv_1/bias/Regularizer/Square/ReadVariableOp�
-Conv_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp6conv_1_bias_regularizer_square_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_1/bias/Regularizer/Square/ReadVariableOp�
Conv_1/bias/Regularizer/SquareSquare5Conv_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_1/bias/Regularizer/Square�
Conv_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_1/bias/Regularizer/Const�
Conv_1/bias/Regularizer/SumSum"Conv_1/bias/Regularizer/Square:y:0&Conv_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/Sum�
Conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_1/bias/Regularizer/mul/x�
Conv_1/bias/Regularizer/mulMul&Conv_1/bias/Regularizer/mul/x:output:0$Conv_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_1/bias/Regularizer/mul�
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
�
c
E__inference_dropout_1_layer_call_and_return_conditional_losses_984128

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:���������&@2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:���������&@2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������&@:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�$
�
B__inference_Conv_3_layer_call_and_return_conditional_losses_985516

inputsA
+conv1d_expanddims_1_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�-Conv_3/bias/Regularizer/Square/ReadVariableOp�/Conv_3/kernel/Regularizer/Square/ReadVariableOp�"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
���������2
conv1d/ExpandDims/dim�
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:���������&@2
conv1d/ExpandDims�
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
conv1d/ExpandDims_1/dim�
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:@ 2
conv1d/ExpandDims_1�
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:���������  *
paddingVALID*
strides
2
conv1d�
conv1d/SqueezeSqueezeconv1d:output:0*
T0*+
_output_shapes
:���������  *
squeeze_dims

���������2
conv1d/Squeeze�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:���������  2	
BiasAdd\
TanhTanhBiasAdd:output:0*
T0*+
_output_shapes
:���������  2
Tanh�
/Conv_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:@ *
dtype021
/Conv_3/kernel/Regularizer/Square/ReadVariableOp�
 Conv_3/kernel/Regularizer/SquareSquare7Conv_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*"
_output_shapes
:@ 2"
 Conv_3/kernel/Regularizer/Square�
Conv_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv_3/kernel/Regularizer/Const�
Conv_3/kernel/Regularizer/SumSum$Conv_3/kernel/Regularizer/Square:y:0(Conv_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/Sum�
Conv_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2!
Conv_3/kernel/Regularizer/mul/x�
Conv_3/kernel/Regularizer/mulMul(Conv_3/kernel/Regularizer/mul/x:output:0&Conv_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/kernel/Regularizer/mul�
-Conv_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02/
-Conv_3/bias/Regularizer/Square/ReadVariableOp�
Conv_3/bias/Regularizer/SquareSquare5Conv_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
: 2 
Conv_3/bias/Regularizer/Square�
Conv_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Conv_3/bias/Regularizer/Const�
Conv_3/bias/Regularizer/SumSum"Conv_3/bias/Regularizer/Square:y:0&Conv_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/Sum�
Conv_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2
Conv_3/bias/Regularizer/mul/x�
Conv_3/bias/Regularizer/mulMul&Conv_3/bias/Regularizer/mul/x:output:0$Conv_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Conv_3/bias/Regularizer/mul�
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp.^Conv_3/bias/Regularizer/Square/ReadVariableOp0^Conv_3/kernel/Regularizer/Square/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*+
_output_shapes
:���������  2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������&@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2^
-Conv_3/bias/Regularizer/Square/ReadVariableOp-Conv_3/bias/Regularizer/Square/ReadVariableOp2b
/Conv_3/kernel/Regularizer/Square/ReadVariableOp/Conv_3/kernel/Regularizer/Square/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:S O
+
_output_shapes
:���������&@
 
_user_specified_nameinputs
�
�
__inference_loss_fn_7_985749E
7dense_1_bias_regularizer_square_readvariableop_resource:@
identity��.Dense_1/bias/Regularizer/Square/ReadVariableOp�
.Dense_1/bias/Regularizer/Square/ReadVariableOpReadVariableOp7dense_1_bias_regularizer_square_readvariableop_resource*
_output_shapes
:@*
dtype020
.Dense_1/bias/Regularizer/Square/ReadVariableOp�
Dense_1/bias/Regularizer/SquareSquare6Dense_1/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2!
Dense_1/bias/Regularizer/Square�
Dense_1/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
Dense_1/bias/Regularizer/Const�
Dense_1/bias/Regularizer/SumSum#Dense_1/bias/Regularizer/Square:y:0'Dense_1/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/Sum�
Dense_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2 
Dense_1/bias/Regularizer/mul/x�
Dense_1/bias/Regularizer/mulMul'Dense_1/bias/Regularizer/mul/x:output:0%Dense_1/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
Dense_1/bias/Regularizer/mul�
IdentityIdentity Dense_1/bias/Regularizer/mul:z:0/^Dense_1/bias/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.Dense_1/bias/Regularizer/Square/ReadVariableOp.Dense_1/bias/Regularizer/Square/ReadVariableOp"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
Z
Convolutional_inputsB
&serving_default_Convolutional_inputs:0����������;
Dense_30
StatefulPartitionedCall:0���������tensorflow/serving/predict:��
�{
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
�__call__
+�&call_and_return_all_conditional_losses
�_default_save_signature"�v
_tf_keras_network�v{"name": "FCN_regressor", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "class_name": "Functional", "config": {"name": "FCN_regressor", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "Convolutional_inputs"}, "name": "Convolutional_inputs", "inbound_nodes": []}, {"class_name": "Conv1D", "config": {"name": "Conv_1", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [25]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_1", "inbound_nodes": [[["Convolutional_inputs", 0, 0, {}]]]}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d", "inbound_nodes": [[["Conv_1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["max_pooling1d", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "Conv_2", "trainable": true, "dtype": "float32", "filters": 64, "kernel_size": {"class_name": "__tuple__", "items": [13]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_2", "inbound_nodes": [[["dropout", 0, 0, {}]]]}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_1", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_1", "inbound_nodes": [[["Conv_2", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout_1", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout_1", "inbound_nodes": [[["max_pooling1d_1", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "Conv_3", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_3", "inbound_nodes": [[["dropout_1", 0, 0, {}]]]}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_2", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_2", "inbound_nodes": [[["Conv_3", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout_2", "inbound_nodes": [[["max_pooling1d_2", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["dropout_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "Dense_1", "trainable": true, "dtype": "float32", "units": 64, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_1", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "Dense_2", "trainable": true, "dtype": "float32", "units": 16, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_2", "inbound_nodes": [[["Dense_1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "Dense_3", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_3", "inbound_nodes": [[["Dense_2", 0, 0, {}]]]}], "input_layers": [["Convolutional_inputs", 0, 0]], "output_layers": [["Dense_3", 0, 0]]}, "shared_object_id": 36, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 200, 8]}, "is_graph_network": true, "save_spec": {"class_name": "TypeSpec", "type_spec": "tf.TensorSpec", "serialized": [{"class_name": "TensorShape", "items": [null, 200, 8]}, "float32", "Convolutional_inputs"]}, "keras_version": "2.5.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "FCN_regressor", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "Convolutional_inputs"}, "name": "Convolutional_inputs", "inbound_nodes": [], "shared_object_id": 0}, {"class_name": "Conv1D", "config": {"name": "Conv_1", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [25]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 1}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 2}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 3}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 4}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_1", "inbound_nodes": [[["Convolutional_inputs", 0, 0, {}]]], "shared_object_id": 5}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d", "inbound_nodes": [[["Conv_1", 0, 0, {}]]], "shared_object_id": 6}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["max_pooling1d", 0, 0, {}]]], "shared_object_id": 7}, {"class_name": "Conv1D", "config": {"name": "Conv_2", "trainable": true, "dtype": "float32", "filters": 64, "kernel_size": {"class_name": "__tuple__", "items": [13]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 8}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 9}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 10}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 11}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_2", "inbound_nodes": [[["dropout", 0, 0, {}]]], "shared_object_id": 12}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_1", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_1", "inbound_nodes": [[["Conv_2", 0, 0, {}]]], "shared_object_id": 13}, {"class_name": "Dropout", "config": {"name": "dropout_1", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout_1", "inbound_nodes": [[["max_pooling1d_1", 0, 0, {}]]], "shared_object_id": 14}, {"class_name": "Conv1D", "config": {"name": "Conv_3", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 15}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 16}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 17}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 18}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Conv_3", "inbound_nodes": [[["dropout_1", 0, 0, {}]]], "shared_object_id": 19}, {"class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_2", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "name": "max_pooling1d_2", "inbound_nodes": [[["Conv_3", 0, 0, {}]]], "shared_object_id": 20}, {"class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout_2", "inbound_nodes": [[["max_pooling1d_2", 0, 0, {}]]], "shared_object_id": 21}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["dropout_2", 0, 0, {}]]], "shared_object_id": 22}, {"class_name": "Dense", "config": {"name": "Dense_1", "trainable": true, "dtype": "float32", "units": 64, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 23}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 24}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 25}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 26}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_1", "inbound_nodes": [[["flatten", 0, 0, {}]]], "shared_object_id": 27}, {"class_name": "Dense", "config": {"name": "Dense_2", "trainable": true, "dtype": "float32", "units": 16, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 28}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 29}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 30}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 31}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_2", "inbound_nodes": [[["Dense_1", 0, 0, {}]]], "shared_object_id": 32}, {"class_name": "Dense", "config": {"name": "Dense_3", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 33}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 34}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "Dense_3", "inbound_nodes": [[["Dense_2", 0, 0, {}]]], "shared_object_id": 35}], "input_layers": [["Convolutional_inputs", 0, 0]], "output_layers": [["Dense_3", 0, 0]]}}, "training_config": {"loss": "mse", "metrics": [[{"class_name": "MeanMetricWrapper", "config": {"name": "mae", "dtype": "float32", "fn": "mean_absolute_error"}, "shared_object_id": 38}]], "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": {"class_name": "CosineDecay", "config": {"initial_learning_rate": 0.001, "decay_steps": 6250, "alpha": 0.05, "name": null}, "shared_object_id": 39}, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
�"�
_tf_keras_input_layer�{"class_name": "InputLayer", "name": "Convolutional_inputs", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 200, 8]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "Convolutional_inputs"}}
�

kernel
bias
	variables
regularization_losses
trainable_variables
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�

_tf_keras_layer�
{"name": "Conv_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Conv1D", "config": {"name": "Conv_1", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [25]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 1}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 2}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 3}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 4}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["Convolutional_inputs", 0, 0, {}]]], "shared_object_id": 5, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 8}}, "shared_object_id": 40}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 200, 8]}}
�
	variables
regularization_losses
trainable_variables
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "max_pooling1d", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "MaxPooling1D", "config": {"name": "max_pooling1d", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "inbound_nodes": [[["Conv_1", 0, 0, {}]]], "shared_object_id": 6, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}, "shared_object_id": 41}}
�
	variables
 regularization_losses
!trainable_variables
"	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "dropout", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "inbound_nodes": [[["max_pooling1d", 0, 0, {}]]], "shared_object_id": 7}
�

#kernel
$bias
%	variables
&regularization_losses
'trainable_variables
(	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�

_tf_keras_layer�
{"name": "Conv_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Conv1D", "config": {"name": "Conv_2", "trainable": true, "dtype": "float32", "filters": 64, "kernel_size": {"class_name": "__tuple__", "items": [13]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 8}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 9}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 10}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 11}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["dropout", 0, 0, {}]]], "shared_object_id": 12, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 32}}, "shared_object_id": 42}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 88, 32]}}
�
)	variables
*regularization_losses
+trainable_variables
,	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "max_pooling1d_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_1", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "inbound_nodes": [[["Conv_2", 0, 0, {}]]], "shared_object_id": 13, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}, "shared_object_id": 43}}
�
-	variables
.regularization_losses
/trainable_variables
0	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "dropout_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dropout", "config": {"name": "dropout_1", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "inbound_nodes": [[["max_pooling1d_1", 0, 0, {}]]], "shared_object_id": 14}
�

1kernel
2bias
3	variables
4regularization_losses
5trainable_variables
6	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�

_tf_keras_layer�
{"name": "Conv_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Conv1D", "config": {"name": "Conv_3", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 15}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 16}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 17}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 18}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["dropout_1", 0, 0, {}]]], "shared_object_id": 19, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 64}}, "shared_object_id": 44}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 38, 64]}}
�
7	variables
8regularization_losses
9trainable_variables
:	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "max_pooling1d_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "MaxPooling1D", "config": {"name": "max_pooling1d_2", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "valid", "data_format": "channels_last"}, "inbound_nodes": [[["Conv_3", 0, 0, {}]]], "shared_object_id": 20, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}, "shared_object_id": 45}}
�
;	variables
<regularization_losses
=trainable_variables
>	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "dropout_2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "inbound_nodes": [[["max_pooling1d_2", 0, 0, {}]]], "shared_object_id": 21}
�
?	variables
@regularization_losses
Atrainable_variables
B	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "flatten", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "inbound_nodes": [[["dropout_2", 0, 0, {}]]], "shared_object_id": 22, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}, "shared_object_id": 46}}
�


Ckernel
Dbias
E	variables
Fregularization_losses
Gtrainable_variables
H	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "Dense_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dense", "config": {"name": "Dense_1", "trainable": true, "dtype": "float32", "units": 64, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 23}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 24}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 25}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 26}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["flatten", 0, 0, {}]]], "shared_object_id": 27, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 512}}, "shared_object_id": 47}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 512]}}
�


Ikernel
Jbias
K	variables
Lregularization_losses
Mtrainable_variables
N	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "Dense_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dense", "config": {"name": "Dense_2", "trainable": true, "dtype": "float32", "units": 16, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 28}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 29}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 30}, "bias_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}, "shared_object_id": 31}, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["Dense_1", 0, 0, {}]]], "shared_object_id": 32, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 64}}, "shared_object_id": 48}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 64]}}
�

Okernel
Pbias
Q	variables
Rregularization_losses
Strainable_variables
T	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"name": "Dense_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "class_name": "Dense", "config": {"name": "Dense_3", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}, "shared_object_id": 33}, "bias_initializer": {"class_name": "Zeros", "config": {}, "shared_object_id": 34}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "inbound_nodes": [[["Dense_2", 0, 0, {}]]], "shared_object_id": 35, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 16}}, "shared_object_id": 49}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 16]}}
�
Uiter

Vbeta_1

Wbeta_2
	Xdecaym�m�#m�$m�1m�2m�Cm�Dm�Im�Jm�Om�Pm�v�v�#v�$v�1v�2v�Cv�Dv�Iv�Jv�Ov�Pv�"
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
p
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9"
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
�
Ylayer_regularization_losses
	variables
Zmetrics

[layers
regularization_losses
trainable_variables
\non_trainable_variables
]layer_metrics
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
#:! 2Conv_1/kernel
: 2Conv_1/bias
.
0
1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
�
^layer_regularization_losses
	variables
_metrics

`layers
regularization_losses
trainable_variables
anon_trainable_variables
blayer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
clayer_regularization_losses
	variables
dmetrics

elayers
regularization_losses
trainable_variables
fnon_trainable_variables
glayer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
hlayer_regularization_losses
	variables
imetrics

jlayers
 regularization_losses
!trainable_variables
knon_trainable_variables
llayer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
#:! @2Conv_2/kernel
:@2Conv_2/bias
.
#0
$1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
�
mlayer_regularization_losses
%	variables
nmetrics

olayers
&regularization_losses
'trainable_variables
pnon_trainable_variables
qlayer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
rlayer_regularization_losses
)	variables
smetrics

tlayers
*regularization_losses
+trainable_variables
unon_trainable_variables
vlayer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
wlayer_regularization_losses
-	variables
xmetrics

ylayers
.regularization_losses
/trainable_variables
znon_trainable_variables
{layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
#:!@ 2Conv_3/kernel
: 2Conv_3/bias
.
10
21"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
.
10
21"
trackable_list_wrapper
�
|layer_regularization_losses
3	variables
}metrics

~layers
4regularization_losses
5trainable_variables
non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
 �layer_regularization_losses
7	variables
�metrics
�layers
8regularization_losses
9trainable_variables
�non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
 �layer_regularization_losses
;	variables
�metrics
�layers
<regularization_losses
=trainable_variables
�non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
 �layer_regularization_losses
?	variables
�metrics
�layers
@regularization_losses
Atrainable_variables
�non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
!:	�@2Dense_1/kernel
:@2Dense_1/bias
.
C0
D1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
.
C0
D1"
trackable_list_wrapper
�
 �layer_regularization_losses
E	variables
�metrics
�layers
Fregularization_losses
Gtrainable_variables
�non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 :@2Dense_2/kernel
:2Dense_2/bias
.
I0
J1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
.
I0
J1"
trackable_list_wrapper
�
 �layer_regularization_losses
K	variables
�metrics
�layers
Lregularization_losses
Mtrainable_variables
�non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 :2Dense_3/kernel
:2Dense_3/bias
.
O0
P1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
O0
P1"
trackable_list_wrapper
�
 �layer_regularization_losses
Q	variables
�metrics
�layers
Rregularization_losses
Strainable_variables
�non_trainable_variables
�layer_metrics
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
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
�0
�1"
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
�0
�1"
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
�0
�1"
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
�0
�1"
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
�0
�1"
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
�

�total

�count
�	variables
�	keras_api"�
_tf_keras_metric�{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}, "shared_object_id": 50}
�

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"�
_tf_keras_metric�{"class_name": "MeanMetricWrapper", "name": "mae", "dtype": "float32", "config": {"name": "mae", "dtype": "float32", "fn": "mean_absolute_error"}, "shared_object_id": 38}
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
(:& 2Adam/Conv_1/kernel/m
: 2Adam/Conv_1/bias/m
(:& @2Adam/Conv_2/kernel/m
:@2Adam/Conv_2/bias/m
(:&@ 2Adam/Conv_3/kernel/m
: 2Adam/Conv_3/bias/m
&:$	�@2Adam/Dense_1/kernel/m
:@2Adam/Dense_1/bias/m
%:#@2Adam/Dense_2/kernel/m
:2Adam/Dense_2/bias/m
%:#2Adam/Dense_3/kernel/m
:2Adam/Dense_3/bias/m
(:& 2Adam/Conv_1/kernel/v
: 2Adam/Conv_1/bias/v
(:& @2Adam/Conv_2/kernel/v
:@2Adam/Conv_2/bias/v
(:&@ 2Adam/Conv_3/kernel/v
: 2Adam/Conv_3/bias/v
&:$	�@2Adam/Dense_1/kernel/v
:@2Adam/Dense_1/bias/v
%:#@2Adam/Dense_2/kernel/v
:2Adam/Dense_2/bias/v
%:#2Adam/Dense_3/kernel/v
:2Adam/Dense_3/bias/v
�2�
.__inference_FCN_regressor_layer_call_fn_984342
.__inference_FCN_regressor_layer_call_fn_984991
.__inference_FCN_regressor_layer_call_fn_985020
.__inference_FCN_regressor_layer_call_fn_984665�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_985157
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_985315
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984766
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984867�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
!__inference__wrapped_model_983994�
���
FullArgSpec
args� 
varargsjargs
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *8�5
3�0
Convolutional_inputs����������
�2�
'__inference_Conv_1_layer_call_fn_985336�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_Conv_1_layer_call_and_return_conditional_losses_985364�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
.__inference_max_pooling1d_layer_call_fn_984009�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *3�0
.�+'���������������������������
�2�
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_984003�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *3�0
.�+'���������������������������
�2�
(__inference_dropout_layer_call_fn_985369
(__inference_dropout_layer_call_fn_985374�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
C__inference_dropout_layer_call_and_return_conditional_losses_985379
C__inference_dropout_layer_call_and_return_conditional_losses_985391�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
'__inference_Conv_2_layer_call_fn_985412�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_Conv_2_layer_call_and_return_conditional_losses_985440�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
0__inference_max_pooling1d_1_layer_call_fn_984024�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *3�0
.�+'���������������������������
�2�
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_984018�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *3�0
.�+'���������������������������
�2�
*__inference_dropout_1_layer_call_fn_985445
*__inference_dropout_1_layer_call_fn_985450�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
E__inference_dropout_1_layer_call_and_return_conditional_losses_985455
E__inference_dropout_1_layer_call_and_return_conditional_losses_985467�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
'__inference_Conv_3_layer_call_fn_985488�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_Conv_3_layer_call_and_return_conditional_losses_985516�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
0__inference_max_pooling1d_2_layer_call_fn_984039�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *3�0
.�+'���������������������������
�2�
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_984033�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *3�0
.�+'���������������������������
�2�
*__inference_dropout_2_layer_call_fn_985521
*__inference_dropout_2_layer_call_fn_985526�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
E__inference_dropout_2_layer_call_and_return_conditional_losses_985531
E__inference_dropout_2_layer_call_and_return_conditional_losses_985543�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
(__inference_flatten_layer_call_fn_985548�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_flatten_layer_call_and_return_conditional_losses_985554�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_Dense_1_layer_call_fn_985575�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_Dense_1_layer_call_and_return_conditional_losses_985598�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_Dense_2_layer_call_fn_985619�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_Dense_2_layer_call_and_return_conditional_losses_985642�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
(__inference_Dense_3_layer_call_fn_985651�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_Dense_3_layer_call_and_return_conditional_losses_985661�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
__inference_loss_fn_0_985672�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_1_985683�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_2_985694�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_3_985705�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_4_985716�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_5_985727�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_6_985738�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_7_985749�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_8_985760�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_9_985771�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�B�
$__inference_signature_wrapper_984962Convolutional_inputs"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 �
B__inference_Conv_1_layer_call_and_return_conditional_losses_985364f4�1
*�'
%�"
inputs����������
� "*�'
 �
0���������� 
� �
'__inference_Conv_1_layer_call_fn_985336Y4�1
*�'
%�"
inputs����������
� "����������� �
B__inference_Conv_2_layer_call_and_return_conditional_losses_985440d#$3�0
)�&
$�!
inputs���������X 
� ")�&
�
0���������L@
� �
'__inference_Conv_2_layer_call_fn_985412W#$3�0
)�&
$�!
inputs���������X 
� "����������L@�
B__inference_Conv_3_layer_call_and_return_conditional_losses_985516d123�0
)�&
$�!
inputs���������&@
� ")�&
�
0���������  
� �
'__inference_Conv_3_layer_call_fn_985488W123�0
)�&
$�!
inputs���������&@
� "����������  �
C__inference_Dense_1_layer_call_and_return_conditional_losses_985598]CD0�-
&�#
!�
inputs����������
� "%�"
�
0���������@
� |
(__inference_Dense_1_layer_call_fn_985575PCD0�-
&�#
!�
inputs����������
� "����������@�
C__inference_Dense_2_layer_call_and_return_conditional_losses_985642\IJ/�,
%�"
 �
inputs���������@
� "%�"
�
0���������
� {
(__inference_Dense_2_layer_call_fn_985619OIJ/�,
%�"
 �
inputs���������@
� "�����������
C__inference_Dense_3_layer_call_and_return_conditional_losses_985661\OP/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� {
(__inference_Dense_3_layer_call_fn_985651OOP/�,
%�"
 �
inputs���������
� "�����������
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984766�#$12CDIJOPJ�G
@�=
3�0
Convolutional_inputs����������
p 

 
� "%�"
�
0���������
� �
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_984867�#$12CDIJOPJ�G
@�=
3�0
Convolutional_inputs����������
p

 
� "%�"
�
0���������
� �
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_985157s#$12CDIJOP<�9
2�/
%�"
inputs����������
p 

 
� "%�"
�
0���������
� �
I__inference_FCN_regressor_layer_call_and_return_conditional_losses_985315s#$12CDIJOP<�9
2�/
%�"
inputs����������
p

 
� "%�"
�
0���������
� �
.__inference_FCN_regressor_layer_call_fn_984342t#$12CDIJOPJ�G
@�=
3�0
Convolutional_inputs����������
p 

 
� "�����������
.__inference_FCN_regressor_layer_call_fn_984665t#$12CDIJOPJ�G
@�=
3�0
Convolutional_inputs����������
p

 
� "�����������
.__inference_FCN_regressor_layer_call_fn_984991f#$12CDIJOP<�9
2�/
%�"
inputs����������
p 

 
� "�����������
.__inference_FCN_regressor_layer_call_fn_985020f#$12CDIJOP<�9
2�/
%�"
inputs����������
p

 
� "�����������
!__inference__wrapped_model_983994�#$12CDIJOPB�?
8�5
3�0
Convolutional_inputs����������
� "1�.
,
Dense_3!�
Dense_3����������
E__inference_dropout_1_layer_call_and_return_conditional_losses_985455d7�4
-�*
$�!
inputs���������&@
p 
� ")�&
�
0���������&@
� �
E__inference_dropout_1_layer_call_and_return_conditional_losses_985467d7�4
-�*
$�!
inputs���������&@
p
� ")�&
�
0���������&@
� �
*__inference_dropout_1_layer_call_fn_985445W7�4
-�*
$�!
inputs���������&@
p 
� "����������&@�
*__inference_dropout_1_layer_call_fn_985450W7�4
-�*
$�!
inputs���������&@
p
� "����������&@�
E__inference_dropout_2_layer_call_and_return_conditional_losses_985531d7�4
-�*
$�!
inputs��������� 
p 
� ")�&
�
0��������� 
� �
E__inference_dropout_2_layer_call_and_return_conditional_losses_985543d7�4
-�*
$�!
inputs��������� 
p
� ")�&
�
0��������� 
� �
*__inference_dropout_2_layer_call_fn_985521W7�4
-�*
$�!
inputs��������� 
p 
� "���������� �
*__inference_dropout_2_layer_call_fn_985526W7�4
-�*
$�!
inputs��������� 
p
� "���������� �
C__inference_dropout_layer_call_and_return_conditional_losses_985379d7�4
-�*
$�!
inputs���������X 
p 
� ")�&
�
0���������X 
� �
C__inference_dropout_layer_call_and_return_conditional_losses_985391d7�4
-�*
$�!
inputs���������X 
p
� ")�&
�
0���������X 
� �
(__inference_dropout_layer_call_fn_985369W7�4
-�*
$�!
inputs���������X 
p 
� "����������X �
(__inference_dropout_layer_call_fn_985374W7�4
-�*
$�!
inputs���������X 
p
� "����������X �
C__inference_flatten_layer_call_and_return_conditional_losses_985554]3�0
)�&
$�!
inputs��������� 
� "&�#
�
0����������
� |
(__inference_flatten_layer_call_fn_985548P3�0
)�&
$�!
inputs��������� 
� "�����������;
__inference_loss_fn_0_985672�

� 
� "� ;
__inference_loss_fn_1_985683�

� 
� "� ;
__inference_loss_fn_2_985694#�

� 
� "� ;
__inference_loss_fn_3_985705$�

� 
� "� ;
__inference_loss_fn_4_9857161�

� 
� "� ;
__inference_loss_fn_5_9857272�

� 
� "� ;
__inference_loss_fn_6_985738C�

� 
� "� ;
__inference_loss_fn_7_985749D�

� 
� "� ;
__inference_loss_fn_8_985760I�

� 
� "� ;
__inference_loss_fn_9_985771J�

� 
� "� �
K__inference_max_pooling1d_1_layer_call_and_return_conditional_losses_984018�E�B
;�8
6�3
inputs'���������������������������
� ";�8
1�.
0'���������������������������
� �
0__inference_max_pooling1d_1_layer_call_fn_984024wE�B
;�8
6�3
inputs'���������������������������
� ".�+'����������������������������
K__inference_max_pooling1d_2_layer_call_and_return_conditional_losses_984033�E�B
;�8
6�3
inputs'���������������������������
� ";�8
1�.
0'���������������������������
� �
0__inference_max_pooling1d_2_layer_call_fn_984039wE�B
;�8
6�3
inputs'���������������������������
� ".�+'����������������������������
I__inference_max_pooling1d_layer_call_and_return_conditional_losses_984003�E�B
;�8
6�3
inputs'���������������������������
� ";�8
1�.
0'���������������������������
� �
.__inference_max_pooling1d_layer_call_fn_984009wE�B
;�8
6�3
inputs'���������������������������
� ".�+'����������������������������
$__inference_signature_wrapper_984962�#$12CDIJOPZ�W
� 
P�M
K
Convolutional_inputs3�0
Convolutional_inputs����������"1�.
,
Dense_3!�
Dense_3���������