
  public void PoissonDisk2DPloy(Polyline poly, double radi, out List<Point3d> pOut, out int ACountOut)
  {
    ///typically k = 30
    int k = 30;

    double cellSize = radi * Math.Sqrt(3);
    Interval cellInterval = new Interval(-1 * cellSize / 2, cellSize / 2);


    List<Box> cells = new List<Box>();
    List<Point3d> activePoints = new List<Point3d>();
    List<Point3d> allPoints = new List<Point3d>();

    Random ran = new Random();
    int initialPointIndex = 0;

    tryNextPoint:

      double xIni = ran.NextDouble() - 0.5;
    double yIni = ran.NextDouble() - 0.5;
    double zIni = ran.NextDouble() - 0.5;


    Point3d initialPoint = new Point3d(poly[initialPointIndex]);

    activePoints.Add(initialPoint);
    //allPoints.Add(initialPoint);


    while(activePoints.Count > 0)
    {

      if( activePoints.Count > 0)
      {

        int randIndex = (int) Math.Floor(ran.NextDouble() * activePoints.Count);
        int pointBadCount = 0;

        for(int t = 0; t < k; t++)
        {
          Point3d pointTemp = NextGaussianSphericalAnnulus(ran, activePoints[randIndex], radi, 1)[0];

          bool insideout = PointInPolygon(pointTemp, poly.ToList());

          if(insideout == false)
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

    if(allPoints.Count == 0 && initialPointIndex < poly.Count - 2)
    {
      initialPointIndex++;
      goto tryNextPoint;
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
      //double zGRSphere = radiRange * zGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);

      PPoints.Add(new Point3d(xGRSphere + Center.X, yGRSphere + Center.Y, 0));
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

  // Function to check if a point is inside a polygon
  static bool PointInPolygon(Point3d point, List<Point3d> polygon)
  {
    int numVertices = polygon.Count;
    double x = point.X;
    double y = point.Y;
    bool inside = false;

    // Store the first point in the polygon and initialize the second point
    Point3d p1 = polygon[0];
    Point3d p2 = Point3d.Unset;

    // Loop through each edge in the polygon
    for (int i = 1; i <= numVertices; i++)
    {
      // Get the next point in the polygon
      p2 = polygon[i % numVertices];

      // Check if the point is above the minimum y coordinate of the edge
      if (y > Math.Min(p1.Y, p2.Y))
      {
        // Check if the point is below the maximum y coordinate of the edge
        if (y <= Math.Max(p1.Y, p2.Y))
        {
          // Check if the point is to the left of the maximum x coordinate of the edge
          if (x <= Math.Max(p1.X, p2.X))
          {
            // Calculate the x-intersection of the line connecting the point to the edge
            double xIntersection = (y - p1.Y) * (p2.X - p1.X) / (p2.Y - p1.Y) + p1.X;

            // Check if the point is on the same line as the edge or to the left of the x-intersection
            if (p1.X == p2.X || x <= xIntersection)
            {
              // Flip the inside flag
              inside = !inside;
            }
          }
        }
      }

      // Store the current point as the first point for the next iteration
      p1 = p2;
    }

    // Return the value of the inside flag
    return inside;
  }
