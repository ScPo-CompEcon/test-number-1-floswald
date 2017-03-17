

module IntTest
	using HW_int
	using Base.Test 

	@testset "check demand" begin
		@test HW_int.q(1) == 2
		@test HW_int.q(4) == 1
		@test HW_int.q(9) == 2/3
	end


	@testset "check gauss adjustment" begin
		@test HW_int.ba2(1,4) == 1.5
		@test HW_int.ab2(1,4) == 2.5
	end

	@testset "eqm condition for Q2" begin
		@test HW_int.dd(1,0,0) == 0
	end
end 