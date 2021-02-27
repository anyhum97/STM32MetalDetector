#define  FT_DIRECT 1

#define  FT_INVERSE 2

static const float Rcoef[14] =
{
	-1.0000000000000000f,  0.0000000000000000f,  0.7071067811865475f,
	0.9238795325112867f,  0.9807852804032304f,  0.9951847266721969f,
	0.9987954562051724f,  0.9996988186962042f,  0.9999247018391445f,
	0.9999811752826011f,  0.9999952938095761f,  0.9999988234517018f,
	0.9999997058628822f,  0.9999999264657178f
};

static const float Icoef[14] =
{
	0.0000000000000000f, -1.0000000000000000f, -0.7071067811865474f,
	-0.3826834323650897f, -0.1950903220161282f, -0.0980171403295606f,
	-0.0490676743274180f, -0.0245412285229122f, -0.0122715382857199f,
	-0.0061358846491544f, -0.0030679567629659f, -0.0015339801862847f,
	-0.0007669903187427f, -0.0003834951875714f
};

int FFT(float* Rdat, float* Idat, int n, int log, int flag)
{
	if(Rdat == 0 || Idat == 0)
	{
		return 0;
	}

	if(n < 16 || n > 16384)
	{
		return 0;
	}

	if(log < 4 || log > 14)
	{
		return 0;
	}

	if((n&(n-1)) > 0)
	{
		return 0;
	}

	if(flag != FT_DIRECT && flag != FT_INVERSE)
	{
		return 0;
	}

	register int iz, jz, nz, kz, io, ie, in, nn;

	float ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

	nn = n >> 1;

	ie = n;

	for(nz=1; nz<=log; ++nz)
	{
		rw = Rcoef[log - nz];
		iw = Icoef[log - nz];

		if(flag == FT_INVERSE)
		{
			iw = -iw;
		}

		in = ie >> 1;

		ru = 1.0f;
		iu = 0.0f;

		for(jz=0; jz<in; ++jz)
		{
			for(iz=jz; iz<n; iz+=ie)
			{
				io = iz + in;

				rtp = Rdat[iz] + Rdat[io];
				itp = Idat[iz] + Idat[io];
				rtq = Rdat[iz] - Rdat[io];
				itq = Idat[iz] - Idat[io];
				Rdat[io] = rtq * ru - itq * iu;
				Idat[io] = itq * ru + rtq * iu;
				Rdat[iz] = rtp;
				Idat[iz] = itp;
			}

			sr = ru;
			ru = ru * rw - iu * iw;
			iu = iu * rw + sr * iw;
		}

		ie >>= 1;
	}

	for(jz=iz=1; iz<n; iz++)
	{
		if(iz < jz)
		{
			io = iz - 1;
			in = jz - 1;
			rtp = Rdat[in];
			itp = Idat[in];
			Rdat[in] = Rdat[io];
			Idat[in] = Idat[io];
			Rdat[io] = rtp;
			Idat[io] = itp;
		}

		kz = nn;

		while(kz < jz)
		{
			jz   = jz - kz;

			kz >>= 1;
		}

		jz = jz + kz;
	}

	if(flag == FT_DIRECT)
	{
		return 1;
	}

	rw = 1.0f / n;

	for(iz=0; iz<n; iz++)
	{
		Rdat[iz] *= rw;
		Idat[iz] *= rw;
	}

	return 1;
}

