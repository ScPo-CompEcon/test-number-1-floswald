

using Base.Test 

@testset "HWintegration Unit Tests" begin
	@testset "testing components" begin
		@testset "check demand" begin
			@test HWintegration.q(1) == 2
			@test HWintegration.q(4) == 1
			@test HWintegration.q(9) == 2/3
		end


		@testset "check gauss adjustment" begin
			@test HWintegration.ba2(1,4) == 1.5
			@test HWintegration.ab2(1,4) == 2.5
		end

		@testset "eqm condition for Q2" begin
			@test HWintegration.dd(1,0,0) == 0
		end
	end

	@testset "Testing Results" begin
		n = 10
		@test HWintegration.question_1b(n)[:I] ≈ HWintegration.A_SOL atol=1e-6
		@test HWintegration.question_1c(n)[:I] ≈ HWintegration.A_SOL atol=0.6
		@test HWintegration.question_1d(n)[:I] ≈ HWintegration.A_SOL atol=0.4
      
		@test HWintegration.question_2a(n)["E[p]"] > 0
		@test HWintegration.question_2a(n)["Var[p]"] > 0
      
		@test HWintegration.question_2b(n)["E[p]"] > 0
		@test HWintegration.question_2b(n)["Var[p]"] > 0

	end
    
end
