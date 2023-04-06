  public List<Point3d> sphPoissonSampling(Point3d Center, double radi = 1, double intensity = 1)
  {
    int numOfPoint = (int) (4 * Math.PI * radi * radi * intensity);

    List<Point3d> PPoints = new List<Point3d>();

    Random ran = new Random();

    for(int i = 0; i < numOfPoint; i++)
    {
      double xGR = NextGaussian(ran);
      double yGR = NextGaussian(ran);
      double zGR = NextGaussian(ran);
      double xGRSphere = radi * xGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      double yGRSphere = radi * yGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      double zGRSphere = radi * zGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);

      PPoints.Add(new Point3d(xGRSphere + Center.X, yGRSphere + Center.Y, zGRSphere + Center.Z));
    }

    return PPoints;
  }

  public double NextGaussian(Random r, double mu = 0, double sigma = 1)
  {
    var u1 = r.NextDouble();
    var u2 = r.NextDouble();

    var rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) *
      Math.Sin(2.0 * Math.PI * u2);

    var rand_normal = mu + sigma * rand_std_normal;

    return rand_normal;
  }
