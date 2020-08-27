package seamcarving;

import astar.WeightedEdge;
import edu.princeton.cs.algs4.Picture;
import seamcarving.util.PictureUtils;
import astar.AStarSolver;
import astar.AStarGraph;
import java.awt.Point;
import java.util.ArrayList;
import java.util.List;
import java.awt.Color;


public class AStarSeamCarver implements SeamCarver {
    private Picture picture;

    public AStarSeamCarver(Picture picture) {
        if (picture == null) {
            throw new NullPointerException("Picture cannot be null.");
        }
        this.picture = new Picture(picture);
    }

    public Picture picture() {
        return new Picture(picture);
    }

    public void setPicture(Picture picture) {
        this.picture = picture;
    }

    public int width() {
        return picture.width();
    }

    public int height() {
        return picture.height();
    }

    public Color get(int x, int y) {
        return picture.get(x, y);
    }

    // Get the RGB value matrix of the processed picture
    public double energy(int x, int y) {
        if (x < 0 || x > width() - 1 || y < 0 || y > height() - 1) {
            throw new IndexOutOfBoundsException();
        }

        int xLeftInput = x - 1;
        int xRightInput = x + 1;

        if (xLeftInput == -1) {
            xLeftInput = width() - 1;
        }

        if (xRightInput == width()) {
            xRightInput = 0;
        }
        Color xLeft = get(xLeftInput, y);
        Color xRight = get(xRightInput, y);
        // x energy
        double xEnergy = Math.pow(xRight.getRed() - xLeft.getRed(), 2) +
                         Math.pow(xRight.getGreen() - xLeft.getGreen(), 2) +
                         Math.pow(xRight.getBlue() - xLeft.getBlue(), 2);

        int yUpInput = y - 1;
        int yDownInput = y + 1;

        if (yUpInput == -1) {
            yUpInput = height() - 1;
        }

        if (yDownInput == height()) {// left x neighbor
            yDownInput = 0;
        }

        Color yUp = get(x, yUpInput);
        Color yDown = get(x, yDownInput);
        // y energy
        double yEnergy = Math.pow(yDown.getRed() - yUp.getRed(), 2) +
                Math.pow(yDown.getGreen() - yUp.getGreen(), 2) +
                Math.pow(yDown.getBlue() - yUp.getBlue(), 2);

        return Math.sqrt(xEnergy + yEnergy);
    }

    public int[] findHorizontalSeam() {
        return new AStarSeamCarver(pictureFlip(picture)).findVerticalSeam();
    }

    public int[] findVerticalSeam() {
        double[][] energyMatrix = PictureUtils.toEnergyMatrix(new AStarSeamCarver(picture));
        AStarSeamGraph graph = new AStarSeamGraph(energyMatrix);
        Point start = new Point(-1, -1);
        Point end = new Point(-2, -2);

        AStarSolver solver = new AStarSolver(graph, start, end, 10);
        List<Point> solution = solver.solution();

        int[] result = new int[solution.size()-2];
        for (int i = 1; i < solution.size()-1; i++) {
            result[i-1] = solution.get(i).x;

        }
        return result;
    }

    private Picture pictureFlip(Picture p) {
        int newHeight = p.width();
        int newWidth = p.height();
        Picture result = new Picture(newWidth, newHeight);
        for (int i = 0; i < newHeight; i++) {
            for (int j = 0; j < newWidth; j++) {
                result.set(j, i, p.get(i, j));
            }
        }
        return result;
    }

    private class AStarSeamGraph implements AStarGraph<Point> {
        private double[][] energyMatrix;
        private Point dummyStart;
        private Point dummyEnd;

        public AStarSeamGraph(double[][] energyMatrix) {
            this.energyMatrix = energyMatrix;
            dummyStart = new Point(-1, -1);
            dummyEnd = new Point(-2, -2);
        }

        @Override
        public List<WeightedEdge<Point>> neighbors(Point v) {
            ArrayList<WeightedEdge<Point>> neighborList = new ArrayList<>();

            if (v.equals(dummyStart)) {
                for (int i = 0; i < width(); i++) {
                    neighborList.add(new WeightedEdge(v, new Point(i, 0), energyMatrix[i][0]));
                }
            } else if (v.getY() == height() - 1){
                neighborList.add(new WeightedEdge(v, dummyEnd, 0));
            } else {
                int xLocation = v.x;
                int yLocation = v.y;

                // Vertical
                Point mid = new Point(xLocation, yLocation + 1);
                Point left = new Point(xLocation - 1, yLocation + 1);
                Point right = new Point(xLocation + 1, yLocation + 1);

                neighborList.add(new WeightedEdge(v, mid, energyMatrix[mid.x][mid.y]));

                if (width() > 1) {
                    if (left.x < 0 && !(right.x > width() - 1)) {
                        neighborList.add(new WeightedEdge(v, right, energyMatrix[right.x][right.y]));
                    } else if (right.x > width() - 1 && !(left.x < 0)) {
                        neighborList.add(new WeightedEdge(v, left, energyMatrix[left.x][left.y]));
                    } else {
                        neighborList.add(new WeightedEdge(v, right, energyMatrix[right.x][right.y]));
                        neighborList.add(new WeightedEdge(v, left, energyMatrix[left.x][left.y]));
                    }
                }
            }
            return neighborList;
        }

        @Override
        public double estimatedDistanceToGoal(Point s, Point goal) {
            return 0;
        }


    }

}
