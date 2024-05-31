private void RunScript(double radi, Mesh meshIn, double radi_Max_0t1, double radi_Min_0t1, ref object pResult, ref object activeP, ref object A)
  {

    List<Point3d> pointResulXY = new List<Point3d>();
    List<Point3d> pointResulYZ = new List<Point3d>();
    List<Point3d> pointResulXZ = new List<Point3d>();
    int ActivePOutXY;
    int ActivePOutYZ;
    int ActivePOutXZ;


    List<Point3d> PTestXY;
    List<Point3d> PTestYZ;
    List<Point3d> PTestXZ;

    PoissonDisk3D_XY(meshIn, radi, 0, radi_Max_0t1, radi_Min_0t1, out pointResulXY, out ActivePOutXY, out PTestXY);
    PoissonDisk3D_YZ(meshIn, radi, 0, radi_Max_0t1, radi_Min_0t1, out pointResulYZ, out ActivePOutYZ, out PTestYZ);
    PoissonDisk3D_XZ(meshIn, radi, 0, radi_Max_0t1, radi_Min_0t1, out pointResulXZ, out ActivePOutXZ, out PTestXZ);

    ////////////Final Filtering
    ///////////////////////////
    List<Point3d> pointResulFinal = new List<Point3d>(pointResulXY);

    for(int i = 0; i < pointResulYZ.Count ; i++)
    {

      bool distanceInCheckTemp = false;

      MeshPoint mPTemp = meshIn.ClosestMeshPoint(pointResulYZ[i], 0.0);
      Color mColor = meshIn.ColorAt(mPTemp);
      double bColor = Convert.ToDouble(Convert.ToInt32(mColor.B));
      double mult = ((radi_Max_0t1 - radi_Min_0t1) * (bColor / 255)) + radi_Min_0t1;

      for(int j = 0; j < pointResulFinal.Count ; j++)
      {

        if(pointResulYZ[i].DistanceTo(pointResulFinal[j]) < radi * mult)
        {
          distanceInCheckTemp = true;
          break;
        }
      }

      if(distanceInCheckTemp == false)
      {
        pointResulFinal.Add(pointResulYZ[i]);
      }
    }

    for(int i = 0; i < pointResulXZ.Count ; i++)
    {

      bool distanceInCheckTemp = false;

      MeshPoint mPTemp = meshIn.ClosestMeshPoint(pointResulXZ[i], 0.0);
      Color mColor = meshIn.ColorAt(mPTemp);
      double bColor = Convert.ToDouble(Convert.ToInt32(mColor.B));
      double mult = ((radi_Max_0t1 - radi_Min_0t1) * (bColor / 255)) + radi_Min_0t1;

      for(int j = 0; j < pointResulFinal.Count ; j++)
      {

        if(pointResulXZ[i].DistanceTo(pointResulFinal[j]) < radi * mult)
        {
          distanceInCheckTemp = true;
          break;
        }
      }

      if(distanceInCheckTemp == false)
      {
        pointResulFinal.Add(pointResulXZ[i]);
      }
    }

    pResult = pointResulFinal;
    activeP = pointResulXZ;
    A = meshIn;
  }

  // <Custom additional code> 
  public void PoissonDisk3D_XY(Mesh meshInput, double radi, int iterTest, double maxMult, double minMult, out List<Point3d> pOut, out int ACountOut, out List<Point3d> PTTT)
  {

    ///typically k = 30
    int k = 30;

    List<Point3d> activePoints = new List<Point3d>();
    List<Point3d> allPoints = new List<Point3d>();

    Random ran = new Random();

    double xIni = ran.NextDouble() - 0.5;
    double yIni = ran.NextDouble() - 0.5;
    double zIni = ran.NextDouble() - 0.5;

    BoundingBox meshBB = meshInput.GetBoundingBox(Plane.WorldXY);

    Point3d minPoint = meshBB.Min;
    Point3d maxPoint = meshBB.Max;
    Point3d midPoint = new Point3d((minPoint.X + maxPoint.X) / 2, (minPoint.Y + maxPoint.Y) / 2, minPoint.Z);

    Rectangle3d recXY = new Rectangle3d(new Plane(midPoint, Vector3d.ZAxis), minPoint, new Point3d(maxPoint.X, maxPoint.Y, minPoint.Z));


    List<Point3d> ptest = new List<Point3d>();

    ////Random Point on Mesh
    int randIndexMeshV = (int) Math.Floor(ran.NextDouble() * meshInput.Vertices.Count);
    Ray3d rayIni = new Ray3d(meshInput.Vertices.Point3dAt(randIndexMeshV), Vector3d.ZAxis);
    double parmIni = Rhino.Geometry.Intersect.Intersection.MeshRay(meshInput, rayIni);

    activePoints.Add(meshInput.Vertices.Point3dAt(randIndexMeshV));
    allPoints.Add(meshInput.Vertices.Point3dAt(randIndexMeshV));

    int iterTemp = 0;

    while(activePoints.Count > 0 && iterTemp < 5000)
    {
      iterTemp++;


      if( activePoints.Count > 0)
      {

        int randIndex = (int) Math.Floor(ran.NextDouble() * activePoints.Count);
        int pointBadCount = 0;

        for(int t = 0; t < k; t++)
        {
          ////Get Mesh Color for Point Gen
          MeshPoint mPTemp = meshInput.ClosestMeshPoint(activePoints[randIndex], 0.0);
          Color mColor = meshInput.ColorAt(mPTemp);
          double bColor = Convert.ToDouble(Convert.ToInt32(mColor.B));
          double mult = ((maxMult - minMult) * (bColor / 255)) + minMult;

          Point3d pointTemp = NextGaussianSphericalAnnulusXY(ran, activePoints[randIndex], minPoint.Z, radi * mult, 1)[0];

          ////Get Raycasted Points
          List<Point3d> pointTempList = InfiniteRayCast(meshInput, pointTemp, Vector3d.ZAxis);
          ptest.Add(pointTemp);

          for(int q = 0; q < pointTempList.Count; q++)
          {
            ////Get Mesh Color for Comparision
            MeshPoint mPTempComp = meshInput.ClosestMeshPoint(pointTempList[q], 0.0);
            Color mColorComp = meshInput.ColorAt(mPTempComp);
            double bColorComp = Convert.ToDouble(Convert.ToInt32(mColorComp.B));
            double multComp = ((maxMult - minMult) * (bColorComp / 255)) + minMult;

            if((pointTempList[q].X < recXY.Center.X - recXY.Width / 2)
              || (pointTempList[q].X > recXY.Center.X + recXY.Width / 2)
              || (pointTempList[q].Y < recXY.Center.Y - recXY.Height / 2)
              || (pointTempList[q].Y > recXY.Center.Y + recXY.Height / 2))
            {
              pointBadCount++;
            }
            else
            {
              bool distanceInCheckTemp = false;

              for(int j = 0; j < allPoints.Count ; j++)
              {

                if(allPoints[j].DistanceTo(pointTempList[q]) < radi * multComp)
                {
                  distanceInCheckTemp = true;
                  break;
                }
              }

              if(distanceInCheckTemp == false)
              {
                activePoints.Add(pointTempList[q]);
                allPoints.Add(pointTempList[q]);
              }
              else
              {
                pointBadCount++;
              }
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
    PTTT = ptest;
  }

  public void PoissonDisk3D_YZ(Mesh meshInput, double radi, int iterTest, double maxMult, double minMult, out List<Point3d> pOut, out int ACountOut, out List<Point3d> PTTT)
  {

    ///typically k = 30
    int k = 30;

    List<Point3d> activePoints = new List<Point3d>();
    List<Point3d> allPoints = new List<Point3d>();

    Random ran = new Random();

    double xIni = ran.NextDouble() - 0.5;
    double yIni = ran.NextDouble() - 0.5;
    double zIni = ran.NextDouble() - 0.5;

    BoundingBox meshBB = meshInput.GetBoundingBox(Plane.WorldXY);

    Point3d minPoint = meshBB.Min;
    Point3d maxPoint = meshBB.Max;
    Point3d midPoint = new Point3d(minPoint.X, (minPoint.Y + maxPoint.Y) / 2, (minPoint.Z + maxPoint.Z) / 2);

    Rectangle3d recYZ = new Rectangle3d(new Plane(midPoint, Vector3d.YAxis), minPoint, new Point3d(maxPoint.X, minPoint.Y, maxPoint.Z));


    List<Point3d> ptest = new List<Point3d>();

    ////Random Point on Mesh
    int randIndexMeshV = (int) Math.Floor(ran.NextDouble() * meshInput.Vertices.Count);
    Ray3d rayIni = new Ray3d(meshInput.Vertices.Point3dAt(randIndexMeshV), Vector3d.ZAxis);
    double parmIni = Rhino.Geometry.Intersect.Intersection.MeshRay(meshInput, rayIni);

    //activePoints.Add(rayIni.PointAt(parmIni));
    //allPoints.Add(rayIni.PointAt(parmIni));

    activePoints.Add(meshInput.Vertices.Point3dAt(randIndexMeshV));
    allPoints.Add(meshInput.Vertices.Point3dAt(randIndexMeshV));

    int iterTemp = 0;

    while(activePoints.Count > 0 && iterTemp < 5000)
    {
      iterTemp++;


      if( activePoints.Count > 0)
      {

        int randIndex = (int) Math.Floor(ran.NextDouble() * activePoints.Count);
        int pointBadCount = 0;

        for(int t = 0; t < k; t++)
        {
          ////Get Mesh Color for Point Gen
          MeshPoint mPTemp = meshInput.ClosestMeshPoint(activePoints[randIndex], 0.0);
          Color mColor = meshInput.ColorAt(mPTemp);
          double bColor = Convert.ToDouble(Convert.ToInt32(mColor.B));
          double mult = ((maxMult - minMult) * (bColor / 255)) + minMult;

          Point3d pointTemp = NextGaussianSphericalAnnulusYZ(ran, activePoints[randIndex], minPoint.X, radi * mult, 1)[0];

          ////Get Raycasted Points
          List<Point3d> pointTempList = InfiniteRayCast(meshInput, pointTemp, Vector3d.XAxis);
          ptest.Add(pointTemp);

          for(int q = 0; q < pointTempList.Count; q++)
          {
            ////Get Mesh Color for Comparision
            MeshPoint mPTempComp = meshInput.ClosestMeshPoint(pointTempList[q], 0.0);
            Color mColorComp = meshInput.ColorAt(mPTempComp);
            double bColorComp = Convert.ToDouble(Convert.ToInt32(mColorComp.B));
            double multComp = ((maxMult - minMult) * (bColorComp / 255)) + minMult;


            if((pointTempList[q].Y < recYZ.Center.Y - (maxPoint.Y - minPoint.Y) / 2)
              || (pointTempList[q].Y > recYZ.Center.Y + (maxPoint.Y - minPoint.Y) / 2)
              || (pointTempList[q].Z < recYZ.Center.Z - (maxPoint.Z - minPoint.Z) / 2)
              || (pointTempList[q].Z > recYZ.Center.Z + (maxPoint.Z - minPoint.Z) / 2))
            {
              pointBadCount++;
            }
            else
            {
              bool distanceInCheckTemp = false;

              for(int j = 0; j < allPoints.Count ; j++)
              {

                if(allPoints[j].DistanceTo(pointTempList[q]) < radi * multComp)
                {
                  distanceInCheckTemp = true;
                  break;
                }
              }

              if(distanceInCheckTemp == false)
              {
                activePoints.Add(pointTempList[q]);
                allPoints.Add(pointTempList[q]);
              }
              else
              {
                pointBadCount++;
              }
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
    PTTT = ptest;
  }

  public void PoissonDisk3D_XZ(Mesh meshInput, double radi, int iterTest, double maxMult, double minMult, out List<Point3d> pOut, out int ACountOut, out List<Point3d> PTTT)
  {

    ///typically k = 30
    int k = 30;

    List<Point3d> activePoints = new List<Point3d>();
    List<Point3d> allPoints = new List<Point3d>();

    Random ran = new Random();

    double xIni = ran.NextDouble() - 0.5;
    double yIni = ran.NextDouble() - 0.5;
    double zIni = ran.NextDouble() - 0.5;

    BoundingBox meshBB = meshInput.GetBoundingBox(Plane.WorldXY);

    Point3d minPoint = meshBB.Min;
    Point3d maxPoint = meshBB.Max;
    Point3d midPoint = new Point3d((minPoint.X + maxPoint.X) / 2, minPoint.Y, (minPoint.Z + maxPoint.Z) / 2);

    Rectangle3d recXZ = new Rectangle3d(new Plane(midPoint, Vector3d.XAxis), minPoint, new Point3d(minPoint.X, maxPoint.Y, maxPoint.Z));


    List<Point3d> ptest = new List<Point3d>();

    ////Random Point on Mesh
    int randIndexMeshV = (int) Math.Floor(ran.NextDouble() * meshInput.Vertices.Count);
    Ray3d rayIni = new Ray3d(meshInput.Vertices.Point3dAt(randIndexMeshV), Vector3d.ZAxis);
    double parmIni = Rhino.Geometry.Intersect.Intersection.MeshRay(meshInput, rayIni);

    //activePoints.Add(rayIni.PointAt(parmIni));
    //allPoints.Add(rayIni.PointAt(parmIni));

    activePoints.Add(meshInput.Vertices.Point3dAt(randIndexMeshV));
    allPoints.Add(meshInput.Vertices.Point3dAt(randIndexMeshV));

    int iterTemp = 0;

    while(activePoints.Count > 0 && iterTemp < 5000)
    {
      iterTemp++;


      if( activePoints.Count > 0)
      {

        int randIndex = (int) Math.Floor(ran.NextDouble() * activePoints.Count);
        int pointBadCount = 0;

        for(int t = 0; t < k; t++)
        {
          ////Get Mesh Color for Point Gen
          MeshPoint mPTemp = meshInput.ClosestMeshPoint(activePoints[randIndex], 0.0);
          Color mColor = meshInput.ColorAt(mPTemp);
          double bColor = Convert.ToDouble(Convert.ToInt32(mColor.B));
          double mult = ((maxMult - minMult) * (bColor / 255)) + minMult;

          Point3d pointTemp = NextGaussianSphericalAnnulusXZ(ran, activePoints[randIndex], minPoint.Y, radi * mult, 1)[0];

          ////Get Raycasted Points
          List<Point3d> pointTempList = InfiniteRayCast(meshInput, pointTemp, Vector3d.YAxis);
          ptest.Add(pointTemp);

          for(int q = 0; q < pointTempList.Count; q++)
          {
            ////Get Mesh Color for Comparision
            MeshPoint mPTempComp = meshInput.ClosestMeshPoint(pointTempList[q], 0.0);
            Color mColorComp = meshInput.ColorAt(mPTempComp);
            double bColorComp = Convert.ToDouble(Convert.ToInt32(mColorComp.B));
            double multComp = ((maxMult - minMult) * (bColorComp / 255)) + minMult;


            if((pointTempList[q].X < recXZ.Center.X - (maxPoint.X - minPoint.X) / 2)
              || (pointTempList[q].X > recXZ.Center.X + (maxPoint.X - minPoint.X) / 2)
              || (pointTempList[q].Z < recXZ.Center.Z - (maxPoint.Z - minPoint.Z) / 2)
              || (pointTempList[q].Z > recXZ.Center.Z + (maxPoint.Z - minPoint.Z) / 2))
            {
              pointBadCount++;
            }
            else
            {
              bool distanceInCheckTemp = false;

              for(int j = 0; j < allPoints.Count ; j++)
              {

                if(allPoints[j].DistanceTo(pointTempList[q]) < radi * multComp)
                {
                  distanceInCheckTemp = true;
                  break;
                }
              }

              if(distanceInCheckTemp == false)
              {
                activePoints.Add(pointTempList[q]);
                allPoints.Add(pointTempList[q]);
              }
              else
              {
                pointBadCount++;
              }
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
    PTTT = ptest;
  }

  public List<Point3d> NextGaussianSphericalAnnulusXY(Random ran, Point3d Center, double ZBottom, double radi = 1, int k = 30)
  {
    int numOfPoint = k;

    List<Point3d> PPoints = new List<Point3d>();

    for(int i = 0; i < numOfPoint; i++)
    {
      double xGR = NextGaussian(ran);
      double yGR = NextGaussian(ran);
      double zGR = NextGaussian(ran);
      double radiRange = radi + ran.NextDouble() * radi;

      double xGRSphere = radiRange * xGR / Math.Sqrt(xGR * xGR + yGR * yGR);
      double yGRSphere = radiRange * yGR / Math.Sqrt(xGR * xGR + yGR * yGR);
      //double xGRSphere = radiRange * xGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      //double yGRSphere = radiRange * yGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      //double zGRSphere = radiRange * zGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);

      PPoints.Add(new Point3d(xGRSphere + Center.X, yGRSphere + Center.Y, ZBottom));
    }

    return PPoints;
  }


  public List<Point3d> NextGaussianSphericalAnnulusYZ(Random ran, Point3d Center, double XBottom, double radi = 1, int k = 30)
  {
    int numOfPoint = k;

    List<Point3d> PPoints = new List<Point3d>();

    for(int i = 0; i < numOfPoint; i++)
    {
      double xGR = NextGaussian(ran);
      double yGR = NextGaussian(ran);
      double zGR = NextGaussian(ran);
      double radiRange = radi + ran.NextDouble() * radi;

      double yGRSphere = radiRange * yGR / Math.Sqrt(yGR * yGR + zGR * zGR);
      double zGRSphere = radiRange * zGR / Math.Sqrt(yGR * yGR + zGR * zGR);
      //double xGRSphere = radiRange * xGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      //double yGRSphere = radiRange * yGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      //double zGRSphere = radiRange * zGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);

      PPoints.Add(new Point3d(XBottom, yGRSphere + Center.Y, zGRSphere + Center.Z));
    }

    return PPoints;
  }

  public List<Point3d> NextGaussianSphericalAnnulusXZ(Random ran, Point3d Center, double YBottom, double radi = 1, int k = 30)
  {
    int numOfPoint = k;

    List<Point3d> PPoints = new List<Point3d>();

    for(int i = 0; i < numOfPoint; i++)
    {
      double xGR = NextGaussian(ran);
      double yGR = NextGaussian(ran);
      double zGR = NextGaussian(ran);
      double radiRange = radi + ran.NextDouble() * radi;

      double xGRSphere = radiRange * xGR / Math.Sqrt(xGR * xGR + zGR * zGR);
      double zGRSphere = radiRange * zGR / Math.Sqrt(xGR * xGR + zGR * zGR);

      //double xGRSphere = radiRange * xGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      //double yGRSphere = radiRange * yGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);
      //double zGRSphere = radiRange * zGR / Math.Sqrt(xGR * xGR + yGR * yGR + zGR * zGR);

      PPoints.Add(new Point3d(xGRSphere + Center.X, YBottom, zGRSphere + Center.Z));
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


  public bool IsPointInPolygon(Point3d[] polygon, Point3d testPoint, double tolerance)
  {
    bool result = false;
    int j = polygon.Length - 1;
    for (int i = 0; i < polygon.Length; i++)
    {
      if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
      {
        if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
        {
          result = !result;
        }
      }

      if(IsInsideLineStrict(new Line(polygon[i], polygon[j]), testPoint.X, testPoint.Y, tolerance))
      {
        return false;
      }
      j = i;
    }

    Array.Clear(polygon, 0, polygon.Length);
    return result;
  }


  private bool IsInsideLineStrict(Line line, double x, double y, double tol)
  {
    return (x >= line.From.X - tol && x <= line.To.X + tol
      || x >= line.To.X - tol && x <= line.From.X + tol)
      && (y >= line.From.Y - tol && y <= line.To.Y + tol
      || y >= line.To.Y - tol && y <= line.From.Y + tol);
  }


  public List<Point3d> InfiniteRayCast(Mesh meshIn, Point3d pointIn, Vector3d vecIn)
  {
    List<Point3d> castedList = new List<Point3d>();

    Ray3d rayIni = new Ray3d(pointIn, vecIn);
    double parmReal = 0;
    double parmTemp = 0;

    int iterLimit = 0;

    while(parmReal >= 0 && iterLimit < 100000)
    {
      rayIni = new Ray3d(rayIni.PointAt(parmTemp), vecIn);
      parmReal = Rhino.Geometry.Intersect.Intersection.MeshRay(meshIn, rayIni);
      parmTemp = parmReal;
      parmTemp += 0.00001;

      if(rayIni.PointAt(parmReal).IsValid)
      {
        castedList.Add(rayIni.PointAt(parmReal));
      }
      else
      {
        break;
      }

      iterLimit++;
    }

    return castedList;
  }
