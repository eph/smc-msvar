module test_msvar
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_swzm1v1

    use model_swzm1v1_t, only: model

    real(wp) :: parasim(100000,54), lpdf

    real(wp) :: p0(54)
    type(model) :: swz
    character(len=144) :: mufile, varfile

    swz = model()

    call assert_equals(54, swz%npara)
    call assert_equals(183, swz%T)
    call assert_equals(3, swz%nobs)
    call assert_equals(5, swz%p)
    call assert_equals(1, swz%ns_mu)
    call assert_equals(1, swz%ns_var)
    call assert_equals(1, swz%ns)


    p0 = [125.730_wp, 84.578_wp, 54.303_wp,-54.529_wp, -3.103_wp, 20.981_wp,-54.890_wp,  8.429_wp,-95.042_wp,142.540_wp, 13.495_wp, 32.189_wp, 13.586_wp, -1.255_wp, 25.950_wp,-39.483_wp, 14.809_wp, 29.036_wp, 17.699_wp, -6.200_wp,  5.008_wp, -0.107_wp,164.214_wp, 99.716_wp, 38.117_wp, 32.166_wp, 19.109_wp,-22.681_wp, 21.814_wp, -7.998_wp,-35.420_wp, 34.484_wp,-25.540_wp, -7.258_wp, 44.071_wp,  6.020_wp, -8.435_wp, -0.025_wp,-247.662_wp, 82.674_wp, 13.246_wp, 10.592_wp,-89.664_wp, -2.935_wp, 17.063_wp, -3.500_wp, -6.658_wp, 18.117_wp,  4.108_wp, 33.954_wp, 17.469_wp,-12.105_wp, -8.910_wp, -0.125_wp]

    call assert_equals(-249.2239_wp, swz%prior%logpdf(p0), 0.001_wp)
    call assert_equals(-3988.8155_wp, swz%lik(p0), 0.02_wp)


  end subroutine test_swzm1v1



  subroutine test_swzm2v1

    use model_swzm2v1_t, only: model

    real(wp) :: parasim(100000,110), lpdf

    real(wp) :: p0(110)
    type(model) :: swz
    character(len=144) :: mufile, varfile

    swz = model()

    call assert_equals(110, swz%npara)
    call assert_equals(183, swz%T)
    call assert_equals(3, swz%nobs)
    call assert_equals(5, swz%p)
    call assert_equals(2, swz%ns_mu)
    call assert_equals(1, swz%ns_var)
    call assert_equals(2, swz%ns)


    p0 = [125.730_wp, 84.578_wp, 54.303_wp,-54.529_wp, -3.103_wp, 20.981_wp,-54.890_wp,  8.429_wp,-95.042_wp,142.540_wp, 13.495_wp, 32.189_wp, 13.586_wp, -1.255_wp, 25.950_wp,-39.483_wp, 14.809_wp, 29.036_wp, 17.699_wp, -6.200_wp,  5.008_wp, -0.107_wp,164.214_wp, 99.716_wp, 38.117_wp, 32.166_wp, 19.109_wp,-22.681_wp, 21.814_wp, -7.998_wp,-35.420_wp, 34.484_wp,-25.540_wp, -7.258_wp, 44.071_wp,  6.020_wp, -8.435_wp, -0.025_wp,-247.662_wp, 82.674_wp, 13.246_wp, 10.592_wp,-89.664_wp, -2.935_wp, 17.063_wp, -3.500_wp, -6.658_wp, 18.117_wp,  4.108_wp, 33.954_wp, 17.469_wp,-12.105_wp, -8.910_wp, -0.125_wp,125.730_wp, 84.578_wp, 54.303_wp,-54.529_wp, -3.103_wp, 20.981_wp,-54.890_wp,  8.429_wp,-95.042_wp,142.540_wp, 13.495_wp, 32.189_wp, 13.586_wp, -1.255_wp, 25.950_wp,-39.483_wp, 14.809_wp, 29.036_wp, 17.699_wp, -6.200_wp,  5.008_wp, -0.107_wp,164.214_wp, 99.716_wp, 38.117_wp, 32.166_wp, 19.109_wp,-22.681_wp, 21.814_wp, -7.998_wp,-35.420_wp, 34.484_wp,-25.540_wp, -7.258_wp, 44.071_wp,  6.020_wp, -8.435_wp, -0.025_wp,-247.662_wp, 82.674_wp, 13.246_wp, 10.592_wp,-89.664_wp, -2.935_wp, 17.063_wp, -3.500_wp, -6.658_wp, 18.117_wp,  4.108_wp, 33.954_wp, 17.469_wp,-12.105_wp, -8.910_wp, -0.125_wp, 0.4_wp, 0.6_wp]

    call assert_equals((-249.2239_wp-2.5416689684691995_wp)*2.0_wp , swz%prior%logpdf(p0), 0.001_wp)
    call assert_equals(-3988.8155_wp, swz%lik(p0), 0.02_wp)


  end subroutine test_swzm2v1


  subroutine test_swzm1v2

    use model_swzm1v2_t, only: model

    real(wp) :: parasim(100000,59), lpdf

    real(wp) :: p0(59)
    type(model) :: swz
    character(len=144) :: mufile, varfile

    swz = model()

    call assert_equals(59, swz%npara)
    call assert_equals(183, swz%T)
    call assert_equals(3, swz%nobs)
    call assert_equals(5, swz%p)
    call assert_equals(1, swz%ns_mu)
    call assert_equals(2, swz%ns_var)
    call assert_equals(2, swz%ns)

    p0 = [125.730_wp, 84.578_wp, 54.303_wp,-54.529_wp, -3.103_wp, 20.981_wp,-54.890_wp,  8.429_wp,-95.042_wp,142.540_wp, 13.495_wp, 32.189_wp, 13.586_wp, -1.255_wp, 25.950_wp,-39.483_wp, 14.809_wp, 29.036_wp, 17.699_wp, -6.200_wp,  5.008_wp, -0.107_wp,164.214_wp, 99.716_wp, 38.117_wp, 32.166_wp, 19.109_wp,-22.681_wp, 21.814_wp, -7.998_wp,-35.420_wp, 34.484_wp,-25.540_wp, -7.258_wp, 44.071_wp,  6.020_wp, -8.435_wp, -0.025_wp,-247.662_wp, 82.674_wp, 13.246_wp, 10.592_wp,-89.664_wp, -2.935_wp, 17.063_wp, -3.500_wp, -6.658_wp, 18.117_wp,  4.108_wp, 33.954_wp, 17.469_wp,-12.105_wp, -8.910_wp, -0.125_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.4_wp, 0.6_wp]

    call assert_equals((-249.2239_wp-2.5416689684691995_wp)*2.0_wp , swz%prior%logpdf(p0), 0.001_wp)
    call assert_equals(-3988.8155_wp, swz%lik(p0), 0.02_wp)


  end subroutine test_swzm1v2


  subroutine test_rfb

    use model_rfbm1v1_t, only: model

    real(wp) :: p0(54)
    type(model) :: rfb

    rfb = model()

    p0 = [-142.2253601596538_wp, -1.49729668218543_wp, -102.4730988377744_wp, -21.684609763784866_wp, -36.553923326408594_wp, 122.85014572464739_wp, -167.69243008534932_wp, -23.571986788112017_wp, -3.898322694372122_wp, 27.0976539288588_wp, 29.255126596177707_wp, 40.80002484452384_wp, -3.786474083534311_wp, -2.450025466546461_wp, -37.89924090238382_wp, -22.19211415743412_wp, 5.818583583883852_wp, 29.444476870060292_wp, 33.22889944891084_wp, -9.563580689223905_wp, -23.455050114439167_wp, -0.0743053131946391_wp, -16.15327192163147_wp, -53.65545897386337_wp, -6.398652148439028_wp, -3.8037320101941985_wp, -5.328659187735242_wp, 5.923710748485433_wp, 10.418276258302427_wp, -21.276867069842694_wp, 1.5191645929354272_wp, -3.7416195393165808_wp, -12.066394499603732_wp, -6.508970318046397_wp, -2.136927291132539_wp, -5.500462913044265_wp, 3.284946481813498_wp, -0.1694919064208663_wp, 11.140861271165333_wp, -4.52463893172639_wp, 112.44280945817403_wp, 21.2291208098362_wp, 14.54251936535313_wp, -40.99168508919068_wp, -44.152753899454446_wp, -5.154880565754422_wp, 41.36938990053708_wp, -15.27547627190866_wp, -18.559002398092186_wp, 13.308939318093996_wp, 7.48400988676366_wp, -1.16873534594915_wp, -15.829352065747853_wp, -0.0625947722115772_wp]

    call assert_equals(54, rfb%npara)
    call assert_equals(183, rfb%T)
    call assert_equals(3, rfb%nobs)
    call assert_equals(5, rfb%p)
    call assert_equals(1, rfb%ns_mu)
    call assert_equals(1, rfb%ns_var)
    call assert_equals(1, rfb%ns)

    ! likelihood

    ! prior 
    call assert_equals(?, rfb%lik(p0))
    call assert_equals(?, rfb%prior%logpdf(p0))


  end subroutine test_rfb

  subroutine test_rfb_hier

    use model_rfb_hierm1v1_t, only: model

    type(model) :: rfb_hier
    real(wp) :: p0(?)

    p0 = [?]

  end subroutine test_rfb_hier

end module test_msvar
