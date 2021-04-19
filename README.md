# delaunator-java
This code is java port of delaunator.

Delaunator triangulation algorithm is incredibly fast and robust library for point of 2D clouds.

I modified it by referring to the code c++, c# version.

to verify contents, go to original site https://github.com/mapbox/delaunator

## Example

```js
public static void main(String[] args) throws IOException
{
    List<DPoint>list = new ArrayList();
    for (int i = 0; i < 10000; i++)
    {
        double x = Math.random() * 1000;
        double y = Math.random() * 1000;
        DPoint p = new DPoint((int)x, (int)y);
        list.add(p);
    }
    
    Delaunator     del = new Delaunator(list.toArray(new DPoint[]{}));
    List<DTriangle>tri = del.getTriangles();
    {
        for (int j = 0; j < tri.size(); j++)
        {
            DTriangle t = tri.get(j);
            DEdge      edgeA = t.ab;        // edge A of Tri
            DEdge      edgeB = t.bc;        // edge B of Tri
            DEdge      edgeC = t.ca;        // edge C of Tri
            DTriangle[]wingA = a.getWing(); // wing A0, A1 of edge A
            DTriangle[]wingB = b.getWing(); // wing B0, B1 of edge B
            DTriangle[]wingC = c.getWing(); // wing C0, C1 of edge C
        }
    }
}
```
