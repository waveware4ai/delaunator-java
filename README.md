# delaunator-java

<img src="delaunator.example.png" alt="delaunay example" width="600" />

This code is java port of delaunator.

Delaunator triangulation algorithm is incredibly fast and robust library for point of 2D clouds.

I modified it by referring to the code c++, c# version.

For more information, go to original site https://github.com/mapbox/delaunator

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
    
    Delaunator     del = new Delaunator(list);
    List<DTriangle>tri = del.getTriangles();
    for (DTriangle t : tri)
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
```
## Ports to other languages

- [delaunator-rs](https://github.com/mourner/delaunator-rs) (Rust)
- [fogleman/delaunay](https://github.com/fogleman/delaunay) (Go)
- [delaunator-cpp](https://github.com/abellgithub/delaunator-cpp) (C++)
- [delaunator-sharp](https://github.com/nol1fe/delaunator-sharp) (C#)
- [delaunator-ruby](https://github.com/hendrixfan/delaunator-ruby) (Ruby)
- [Delaunator-Python](https://github.com/HakanSeven12/Delaunator-Python) (Python)
- [hx-delaunator](https://github.com/dmitryhryppa/hx-delaunator) (Haxe)
- [ricardomatias/delaunator](https://github.com/ricardomatias/delaunator) (Kotlin)
- [delaunator-java](https://github.com/waveware4ai/delaunator-java) (Java)
