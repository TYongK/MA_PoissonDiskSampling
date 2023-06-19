public void PoissonDiskCube(Box cube, double radi, int iterTest, out List<Point3d> pOut, out int ACountOut)
  {
    ///typically k = 30
    int k = 30;

    double cellSize = radi * Math.Sqrt(3);
    Interval cellInterval = new Interval(-1 * cellSize / 2, cellSize / 2);


    List<Box> cells = new List<Box>();
    List<Point3d> activePoints = new List<Point3d>();
    List<Point3d> allPoints = new List<Point3d>();

    Random ran = new Random();

    double xIni = ran.NextDouble() - 0.5;
    double yIni = ran.NextDouble() - 0.5;
    double zIni = ran.NextDouble() - 0.5;

    Point3d initialPoint = new Point3d(cube.Center.X + (xIni * cube.X.Length), cube.Center.Y + (yIni * cube.Y.Length), cube.Center.Z + (zIni * cube.Z.Length));

    activePoints.Add(initialPoint);
    allPoints.Add(initialPoint);


    int iterTemp = 0;

    while(activePoints.Count > 0)
    {
      iterTemp++;


      if( activePoints.Count > 0)
      {

        int randIndex = (int) Math.Floor(ran.NextDouble() * activePoints.Count);
        int pointBadCount = 0;

        for(int t = 0; t < k; t++)
        {
          Point3d pointTemp = NextGaussianSphericalAnnulus(ran, activePoints[randIndex], radi, 1)[0];


          if((pointTemp.X < cube.Center.X - cube.X.Length / 2)
            || (pointTemp.X > cube.Center.X + cube.X.Length / 2)
            || (pointTemp.Y < cube.Center.Y - cube.Y.Length / 2)
            || (pointTemp.Y > cube.Center.Y + cube.Y.Length / 2)
            || (pointTemp.Z < cube.Center.Z - cube.Z.Length / 2)
            || (pointTemp.Z > cube.Center.Z + cube.Z.Length / 2))
          {
            pointBadCount++;
          }
          else
          {

            bool distanceInCheckTemp = false;

            for(int j = 0; j < allPoints.Count ; j++)
            {
              if(allPoints[j].DistanceTo(pointTemp) < radi)
              {
                distanceInCheckTemp = true;
                break;
              }
            }

            if(distanceInCheckTemp == false)
            {
              activePoints.Add(pointTemp);
              allPoints.Add(pointTemp);
            }
            else
            {
              pointBadCount++;
            }
          }
        }

        if(pointBadCount >= k)
        {
          activePoints.RemoveAt(randIndex);
        }
      }

    }

    pOut = allPoints;
    ACountOut = activePoints.Count;
  }

  public List<Point3d> NextGaussianSphericalAnnulus(Random ran, Point3d Center, double radi = 1, int k = 30)
  {
    int numOfPoint = k;

    List<Point3d> PPoints = new List<Point3d>();

    for(int i = 0; i < numOfPoint; i++)
    {
      double xGR = NextGaussian(ran);
      double yGR = NextGaussian(ran);
      double zGR = NextGaussian(ran);
      double radiRange = radi + ran.NextDouble() * radi;

      double xGRSphere = radiRange * xGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      double yGRSphere = radiRange * yGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      double zGRSphere = radiRange * zGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);

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
