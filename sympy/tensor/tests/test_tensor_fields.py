from sympy import Matrix
#from sympy.utilities.pytest import raises

from sympy.tensor import *
from sympy import sin, symbols, cos

def test_df_varlist():
    x1, x2, x3 = symbols('x1 x2 x3')
    f=x1**2*x2+sin(x2*x3-x2)
        
    #var-list
    var_list=[x1, x2, x3]
    
    print ('test_df_varlist  <============================ actual test code')
    assert df(f, var_list) == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    print ('test_df_varlist_type  <============================ actual test code')
    assert isinstance(df(f, var_list), list)    

    print ('test_df_varlist_l  <============================ actual test code')
    assert df(f, var_list, 'l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    print ('test_df_varlist_l_type  <============================ actual test code')
    assert isinstance(df(f, var_list, 'l'), list)

    print ('test_df_varlist_a  <============================ actual test code')
    assert df(f, var_list, 'a') == list2arraypy([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    print ('test_df_varlist_a_type  <============================ actual test code')
    assert isinstance(df(f, var_list, 'a'), Arraypy)
    
    print ('test_df_varlist_t  <============================ actual test code')
    assert df(f, var_list, 't') == list2tensor([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]) 
    assert isinstance(df(f, var_list, 't'),Tensor)
    assert df(f, var_list, 't').type_pq == (0,1)
    
def test_df_var_tnsr0():
    x1, x2, x3 = symbols('x1 x2 x3')
    f=x1**2*x2+sin(x2*x3-x2)            
                
    #var-tensor 0
    var_tnsr0=Tensor(Arraypy(3),(1))
    var_tnsr0[0]=x1
    var_tnsr0[1]=x2
    var_tnsr0[2]=x3 
    
    print ('test_df_var_tnsr0  <============================ actual test code')
    assert df(f, var_tnsr0) == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    assert isinstance(df(f, var_tnsr0), list)     
        
    print ('test_df_var_tnsr0_l  <============================ actual test code')
    assert df(f, var_tnsr0, 'l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    assert isinstance(df(f, var_tnsr0, 'l'), list)    
    
    print ('test_df_var_tnsr0_a  <============================ actual test code')
    assert df(f, var_tnsr0, 'a') == list2arraypy([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    print ('test_df_var_tnsr0_a_type  <============================ actual test code')
    assert isinstance(df(f, var_tnsr0, 'a'), Arraypy)
    
    print ('test_df_var_tnsr0_t  <============================ actual test code')
    assert df(f, var_tnsr0, 't') == list2tensor([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    print ('test_df_var_tnsr0_t_type  <============================ actual test code')
    assert isinstance(df(f, var_tnsr0, 't'), Tensor)
    assert df(f, var_tnsr0, 't').type_pq == (0,1)  
      
def test_df_var_tnsr1():
        x1, x2, x3 = symbols('x1 x2 x3')
        f=x1**2*x2+sin(x2*x3-x2)            
                    
        #var-tensor 1
        var_tnsr1=Arraypy([1,3,1]).to_tensor(1)
        var_tnsr1[1]=x1
        var_tnsr1[2]=x2
        var_tnsr1[3]=x3 
        
        #res-tensor 1
        res_ar1=Arraypy([1,3,1])
        res_ar1[1]=2*x1*x2
        res_ar1[2]=x1**2 + (x3 - 1)*cos(x2*x3 - x2)
        res_ar1[3]=x2*cos(x2*x3 - x2)
        res_ten1 = res_ar1.to_tensor(-1)
        
        print ('test_df_var_tnsr1  <============================ actual test code')
        assert df(f, var_tnsr1) == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
        assert isinstance(df(f, var_tnsr1), list)        
            
        print ('test_df_var_tnsr1_l  <============================ actual test code')
        assert df(f, var_tnsr1, 'l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
        assert isinstance(df(f, var_tnsr1, 'l'), list)    
        
        
        print ('test_df_var_tnsr1_a  <============================ actual test code')
        assert df(f, var_tnsr1, 'a') == res_ar1
        print ('test_df_var_tnsr1_a_type  <============================ actual test code')
        assert isinstance(df(f, var_tnsr1, 'a'), Arraypy)
        
        print ('test_df_var_tnsr1_t  <============================ actual test code')
        assert df(f, var_tnsr1, 't') == res_ten1
        print ('test_df_var_tnsr1_t_type  <============================ actual test code')
        assert isinstance(df(f, var_tnsr1, 't'), Tensor)
        assert df(f, var_tnsr1, 't').type_pq == (0,1)
        
        
def test_div_var_x_list():
    
    x1, x2, x3 = symbols('x1 x2 x3')
    
    X = [x1*x2**3,x2-cos(x3),x3**3-x1]
    var= [x1, x2, x3]
    g=Matrix([[2,1,0],[1,3,0],[0,0,1]])
    
    ten=Arraypy([2,3,0]).to_tensor((-1,-1))
    ten[0,0]=2
    ten[0,1]=1
    ten[0,2]=0
    ten[1,0]=1
    ten[1,1]=3
    ten[1,2]=0
    ten[2,0]=0
    ten[2,1]=0
    ten[2,2]=1
    
    #g задано tensor, индекс с 1 и var-list
    ten1=Arraypy([2,3,1]).to_tensor((-1,-1))
    ten1[1,1]=2
    ten1[1,2]=1
    ten1[1,3]=0
    ten1[2,1]=1
    ten1[2,2]=3
    ten1[2,3]=0
    ten1[3,1]=0
    ten1[3,2]=0
    ten1[3,3]=1 
    
    print ('test_div_var_x_list  <============================ actual test code')
    assert div(X, var) == x2**3 + 3*x3**2 + 1     
    
    print ('test_div_var_x_list_g  <============================ actual test code')
    assert div(X, var, g) == x2**3 + 3*x3**2 + 1     
    
    print ('test_div_var_x_list_t0  <============================ actual test code')
    assert div(X, var, ten) == x2**3 + 3*x3**2 + 1  
    
    print ('test_div_var_x_list_t1  <============================ actual test code')
    assert div(X, var, ten1) == x2**3 + 3*x3**2 + 1 
    
def test_grad_varlist():
    x1, x2, x3= symbols('x1 x2 x3')
    
    f=x1**2*x2+sin(x2*x3-x2)
    var1=[x1, x2, x3]

    res_ar1=Arraypy([1,3,0])
    res_ar1[0]=2*x1*x2
    res_ar1[1]=x1**2 + (x3 - 1)*cos(x2*x3 - x2)
    res_ar1[2]=x2*cos(x2*x3 - x2)
    res_ten1 = res_ar1.to_tensor(1)
    
    res_ar=Arraypy([1,3,0])
    res_ar[0]=-x1**2/5 + 6*x1*x2/5 - (x3 - 1)*cos(x2*x3 - x2)/5
    res_ar[1]=2*x1**2/5 - 2*x1*x2/5 + cos(x2*x3 - x2)*2*(x3 - 1)/5
    res_ar[2]=x2*cos(x2*x3 - x2)
    res_ten = res_ar.to_tensor(1)    
    
    g=Matrix([[2,1,0],[1,3,0],[0,0,1]])
    
    print ('test_grad_l <============================ actual test code')
    assert grad(f,var1,output_type='l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    assert isinstance(grad(f,var1,output_type='l'), list) 
    
    print ('test_grad_a  <============================ actual test code')
    assert grad(f,var1,output_type='a') == res_ar1
    assert isinstance(grad(f,var1,output_type='t'), Arraypy)    

    print ('test_grad_t  <============================ actual test code')
    assert grad(f,var1,output_type='t') == res_ten1
    assert isinstance(grad(f,var1,output_type='t'), Tensor)    
    assert grad(f,var1,output_type='t').type_pq == (1,0)
    
    
    print ('test_grad_g_l <============================ actual test code')
    assert str(grad(f,var1, g, output_type='l')) == '[-x1**2/5 + 6*x1*x2/5 - (x3 - 1)*cos(x2*x3 - x2)/5, 2*x1**2/5 - 2*x1*x2/5 + 2*(x3 - 1)*cos(x2*x3 - x2)/5, x2*cos(x2*x3 - x2)]'
    assert isinstance(grad(f,var1, g, output_type='l'), list)
    
    print ('test_grad_g_a <============================ actual test code')
    assert grad(f,var1, g, output_type='a') == res_ar 
    assert isinstance(grad(f,var1, g, output_type='a'), Arraypy)   
    
    print ('test_grad_g_t <============================ actual test code')
    assert grad(f,var1, g, output_type='t') == res_ten 
    assert isinstance(grad(f,var1, g, output_type='t'), Tensor)    
    assert grad(f,var1, g, output_type='t').type_pq == (1,0) 
    
def test_grad_gtnsr():
    x1, x2, x3= symbols('x1 x2 x3')  
    f=x1**2*x2+sin(x2*x3-x2)
    var1=[x1, x2, x3]    
    
    k1=Arraypy([1,3,0]).to_tensor(1)
    k1[0]=x1
    k1[1]=x2
    k1[2]=x3   
    
    #g задано tensor, индекс с 1 и var-list
    a=Arraypy([2,3,1])
    b=a.to_tensor((-1,-1))
    b[1,1]=2
    b[1,2]=1
    b[1,3]=0
    b[2,1]=1
    b[2,2]=3
    b[2,3]=0
    b[3,1]=0
    b[3,2]=0
    b[3,3]=1   
    
    res_ar=Arraypy([1,3,1])
    res_ar[1]=-x1**2/5 + 6*x1*x2/5 - (x3 - 1)*cos(x2*x3 - x2)/5
    res_ar[2]=2*x1**2/5 - 2*x1*x2/5 + cos(x2*x3 - x2)*2*(x3 - 1)/5
    res_ar[3]=x2*cos(x2*x3 - x2)
    res_ten = res_ar.to_tensor(1) 
    
    print ('test_grad_g_a <============================ actual test code')
    assert grad(f,var1,b,'a') == res_ar 
    assert isinstance(grad(f,var1,b,'a'), Arraypy)    
    
    print ('test_grad_vartnsr_t <============================ actual test code')
    assert grad(f,var1,b,'t') == res_ten 
    assert isinstance(grad(f,var1,b,'t'), Tensor)    
    assert grad(f,var1,b,'t').type_pq == (1,0)
    
def test_LieXY():
    x1, x2, x3,t,l,a = symbols('x1 x2 x3 t l a')
    
    X=[x1*x2**3,x2-cos(x3),x3**3-x1]
    Y=[x1**3*x2**3,x2*x3-sin(x1*x3),x3**3-x1**2]
    arg=[x1, x2, x3] 
    
    res_ar=Arraypy([1,3,0])
    res_ar[0]=2*x1**3*x2**6 + 3*x1**3*x2**2*(x2 - cos(x3)) - 3*x1*x2**2*(x2*x3 - sin(x1*x3))
    res_ar[1]=-x1*x2**3*x3*cos(x1*x3) - x2*x3 + x3*(x2 - cos(x3)) + (-x1 + x3**3)*(-x1*cos(x1*x3) + x2) - (-x1**2 + x3**3)*sin(x3) + sin(x1*x3)
    res_ar[2]=x1**3*x2**3 - 2*x1**2*x2**3 + 3*x3**2*(-x1 + x3**3) - 3*x3**2*(-x1**2 + x3**3)
    res_ten = res_ar.to_tensor(1)
    
    print ('test_LieXY <============================ actual test code')
    assert LieXY(X,Y,arg,'l') == [2*x1**3*x2**6 + 3*x1**3*x2**2*(x2 - cos(x3)) - 3*x1*x2**2*(x2*x3 - sin(x1*x3)), -x1*x2**3*x3*cos(x1*x3) - x2*x3 + x3*(x2 - cos(x3)) + (-x1 + x3**3)*(-x1*cos(x1*x3) + x2) - (-x1**2 + x3**3)*sin(x3) + sin(x1*x3), x1**3*x2**3 - 2*x1**2*x2**3 + 3*x3**2*(-x1 + x3**3) - 3*x3**2*(-x1**2 + x3**3)]
    assert isinstance(LieXY(X,Y,arg,'l'), list)     
    
    
    print ('test_LieXY <============================ actual test code')
    assert LieXY(X,Y,arg,'a') == res_ar
    assert isinstance(LieXY(X,Y,arg,'a'), Arraypy)     
    
    print ('test_LieXY <============================ actual test code')
    assert LieXY(X,Y,arg,'t') == res_ten 
    assert isinstance(LieXY(X,Y,arg,'t'), Tensor)    
    assert LieXY(X,Y,arg,'t').type_pq == (1,0) 

def test_rot():
    x1, x2, x3= symbols('x1 x2 x3')
    X=[x1*x2**3,x2-cos(x3),x3**3-x1]
    arg=[x1,x2,x3]
    
    #X-tensor (0)
    j=Arraypy(3)
    k0=Tensor(j,(1))
    k0[0]=x1*x2**3
    k0[1]=x2-cos(x3)
    k0[2]=x3**3-x1
    
    #X-tensor (1)
    k1=Arraypy([1,3,1]).to_tensor(1)
    k1[1]=x1*x2**3
    k1[2]=x2-cos(x3)
    k1[3]=x3**3-x1  
    
    #args-ten(0)
    v0=Tensor(j,(1))
    v0[0]=x1
    v0[1]=x2
    v0[2]=x3
    
    #args-tensor (1)
    v1=Arraypy([1,3,1]).to_tensor(1)
    v1[1]=x1
    v1[2]=x2
    v1[3]=x3 
    
    s0=Tensor(j,(1))
    s0[0]=-sin(x3)
    s0[1]=1
    s0[2]=-3*x1*x2**2    
    
    print ('test_rot_l <============================ actual test code')
    assert rot(X,arg) == [-sin(x3), 1, -3*x1*x2**2]
    assert isinstance(rot(X,arg), list)
    
    print ('test_rot_a <============================ actual test code')
    assert rot(X,arg, 'a') == list2arraypy([-sin(x3), 1, -3*x1*x2**2])
    assert isinstance(rot(X,arg, 'a'), Arraypy)   
    
    print ('test_rot_t <============================ actual test code')
    assert rot(X,arg, 't') == s0
    assert isinstance(rot(X,arg, 't'), Tensor)     
    assert rot(X,arg, 't').type_pq == (1,0)
    
    
    print ('test_rot_Xt_l <============================ actual test code')
    assert rot(k0,arg) == s0
    assert isinstance(rot(k0,arg), Tensor)     
    assert rot(X,arg, 't').type_pq == (1,0)
    
    print ('test_rot_Xt_a <============================ actual test code')
    assert rot(k0,arg, 'a') == list2arraypy([-sin(x3), 1, -3*x1*x2**2])
    assert isinstance(rot(k0,arg, 'a'), Arraypy)   
        
    print ('test_rot_t <============================ actual test code')
    assert rot(k0,arg, 't') == s0
    assert isinstance(rot(k0,arg, 't'), Tensor)     
    assert rot(X,arg, 't').type_pq == (1,0)    
    
    
    
def test_Lie_w():
    
    x1,x2,x3= symbols('x1, x2, x3')
        
    X=[x1*x2**3,x2-cos(x3),x3**3-x1]
    
    arr=Arraypy((3,3))
    y = Tensor(arr, (-1, -1))
    y1 = Tensor(arr, (-1, -1))
    
    y[0,1]=x3
    y[0,2]=-x2 
    y[1,0]=-x3
    y[1,2]=x1
    y[2,0]=x2 
    y[2,1]=-x1 
    
    y1[0,1]=x2**3*x3 + x3**3 + x3
    y1[0,2]=-x2**4 - 3*x2*x3**2 - x2 + x3*sin(x3) + cos(x3) 
    y1[1,0]=-x2**3*x3 - x3**3 - x3
    y1[1,2]=-2*x1*x2**3 + 3*x1*x3**2 + x1
    y1[2,0]=x2**4 + 3*x2*x3**2 + x2 - x3*sin(x3) - cos(x3) 
    y1[2,1]=2*x1*x2**3 - 3*x1*x3**2 - x1  
    
    """
    print ('test_Lie_w <============================ actual test code')
    assert Lie_w(y, X, [x1,x2,x3]) == y1 
    assert isinstance(Lie_w(y, X, [x1,x2,x3]), Tensor)    
    assert Lie_w(y, X, [x1,x2,x3]).type_pq == (0,2) 
    """
    
    omega = Tensor(arr, (-1, -1))
    omega[0,1]=x2
    omega[1,0]=-x2
    omega[0,2]=-x1 
    omega[2,0]=x1 
    
    ten = Tensor(arr, (-1, -1))
    ten [0,1]=x2**4 + 2*x2 - cos(x3)
    ten [0,2]=-2*x1*x2**3 - 3*x1*x3**2 + x2*sin(x3) 
    ten [1,0]=-x2**4 - 2*x2 + cos(x3)
    ten [1,2]=-3*x1**2*x2**2
    ten [2,0]=2*x1*x2**3 + 3*x1*x3**2 - x2*sin(x3) 
    ten [2,1]=3*x1**2*x2**2 
    
    """
    print ('test_Lie_w <============================ actual test code')
    assert Lie_w(omega, X, [x1,x2,x3]) == ten  
    assert isinstance(Lie_w(omega, X, [x1,x2,x3]), Tensor)    
    assert Lie_w(omega, X, [x1,x2,x3]).type_pq == (0,2) 
    """  

def test_dw():
    x1, x2, x3 = symbols('x1 x2 x3')
    
    y = Tensor(Arraypy((3,3)), (-1,-1))
    y1=Tensor(Arraypy((3,3,3)), (-1,-1,-1))
    
    y[0,1]=x3
    y[0,2]=-x2 
    y[1,0]=-x3
    y[1,2]=x1
    y[2,0]=x2 
    y[2,1]=-x1
    
    y1[0,1,2]=3
    y1[0,2,1]=-3 
    y1[1,0,2]=-3 
    y1[1,2,0]=3
    y1[2,0,1]=3
    y1[2,1,0]=-3
    
    """
    print ('test_dw <============================ actual test code')
    assert dw(y, [x1,x2,x3]) == y1 
    assert isinstance(dw(y, [x1,x2,x3]), Tensor)    
    assert dw(y, [x1,x2,x3]).type_pq == (0,3) 
    """      
    
