
using ForwardDiff
using JSMDUtils.Autodiff

@testset "Autodiff" verbose=true begin 

    fun1(x) = 3x^2 + 4*cos(x) + x*sqrt(x)
    fun2(x) = 2*sin(x[1]) + 3*x[2]^2
    fun3(x) = [2cos(x[1]) + x[2]^2, sin(x[1])^2 + x[2]^3, x[2]/8]
    
    function fun3!(y, x)
        y[1] = 2cos(x[1]) + x[2]^2
        y[2] = sin(x[1])^2 + x[2]^3
        y[3] = x[2]/8
    end

    atol, rtol = 1e-12, 1e-12 

    for _ in 1:10

        xs = rand() 
        xv = rand(2)
        yv = zeros(3)
        
        # -----------
        # Derivatives

        ad = Autodiff.derivative(fun1, xs)
        af = ForwardDiff.derivative(fun1, xs)
        @test ad ≈ af atol=atol rtol=rtol 

        ad = Autodiff.derivative(x->Autodiff.derivative(fun1, x), xs)
        af = JSMDUtils.Math.D²(fun1, xs)
        @test ad ≈ af atol=atol rtol=rtol


        # -----------
        # Gradients 
        a_resg = zeros(size(xv))
        cfg = Autodiff.AutoGradientConfig(xv)

        f_resg = ForwardDiff.gradient(fun2, xv)

        Autodiff.gradient!(a_resg, fun2, xv, cfg)
        @test a_resg ≈ f_resg atol=atol rtol=rtol

        Autodiff.gradient!(a_resg, fun2, xv)
        @test a_resg ≈ f_resg atol=atol rtol=rtol


        # -----------
        # Jacobians

        a_resj = zeros(length(yv), length(xv))
        cfj = Autodiff.AutoJacobianConfig(xv)
        cfj! = Autodiff.AutoJacobianConfig(yv, xv)

        f_resj = ForwardDiff.jacobian(fun3, xv)
        
        Autodiff.jacobian!(a_resj, fun3, xv, cfj)
        @test a_resj ≈ f_resj atol=atol rtol=rtol

        Autodiff.jacobian!(a_resj, fun3, xv)
        @test a_resj ≈ f_resj atol=atol rtol=rtol 

        Autodiff.jacobian!(a_resj, fun3!, yv, xv, cfj!)
        @test a_resj ≈ f_resj atol=atol rtol=rtol 

        Autodiff.jacobian!(a_resj, fun3!, yv, xv)
        @test a_resj ≈ f_resj atol=atol rtol=rtol


        # -----------
        # Hessian 
        a_resh = zeros(length(xv), length(xv))
        cfh = Autodiff.AutoHessianConfig(xv)

        f_resh = ForwardDiff.hessian(fun2, xv)

        Autodiff.hessian!(a_resh, fun2, xv)
        @test a_resh ≈ f_resh atol=atol rtol=rtol

        Autodiff.hessian!(a_resh, fun2, xv, cfh)
        @test a_resh ≈ f_resh atol=atol rtol=rtol

    end

end;
