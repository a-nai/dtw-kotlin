import android.Manifest
import android.annotation.SuppressLint
import android.content.Context
import android.media.AudioManager
import android.media.MediaPlayer
import android.media.MediaRecorder
import android.os.Bundle
import android.os.Environment
import android.view.Menu
import android.view.MenuItem
import android.widget.Button
import android.widget.TextView
import androidx.appcompat.app.AppCompatActivity
import androidx.core.app.ActivityCompat
import com.example.musict.Arrays.contains
import com.example.musict.DTW.getWarpInfoBetween
import com.example.musict.TypeConversions.doubleArrayToByteArray
import com.google.android.material.snackbar.Snackbar
import java.io.*
import java.math.BigInteger
import java.util.*
import java.util.Arrays
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.abs
import kotlin.math.asin
import kotlin.math.pow
import kotlin.math.sqrt


class LinearWindow @JvmOverloads constructor(
    tsI: TimeSeries,
    tsJ: TimeSeries,
    searchRadius: Int = DEFAULT_RADIUS
) :
    SearchWindow(tsI.size(), tsJ.size()) {
    companion object {
        // CONSTANTS
        private const val DEFAULT_RADIUS = 0
    }

    // CONSTRUCTORS
    init {
        val ijRatio = tsI.size().toDouble() / tsJ.size().toDouble()
        val isIlargest = tsI.size() >= tsJ.size()
        for (i in 0 until tsI.size()) {
            if (isIlargest) {
                val j =
                    Math.min(Math.round(i / ijRatio).toInt(), tsJ.size() - 1)
                super.markVisited(i, j)
            } else {
                val maxJ = Math.round((i + 1) / ijRatio).toInt() - 1
                val minJ = Math.round(i / ijRatio).toInt()
                super.markVisited(i, minJ)
                super.markVisited(i, maxJ)
            } // end if
        } // end for loop
        super.expandWindow(searchRadius)
    } // end Constructor
} // end class FullWindow


class FullWindow(tsI: TimeSeries, tsJ: TimeSeries) :
    SearchWindow(tsI.size(), tsJ.size()) {
    // CONSTRUCTOR
    init {
        for (i in 0 until tsI.size()) {
            super.markVisited(i, minJ())
            super.markVisited(i, maxJ())
        } // end for loop
    } // end CONSTRUCTOR
} // end class FullWindow


internal class MemoryResidentMatrix(// PRIVATE DATA
    private val window: SearchWindow
) : CostMatrix {
    private val cellValues: DoubleArray
    private val colOffsets: IntArray

    // PUBLIC FUNCTIONS
    override fun put(col: Int, row: Int, value: Double) {
        if (row < window.minJforI(col) || row > window.maxJforI(col)) {
            throw InternalError(
                "CostMatrix is filled in a cell (col=" + col + ", row=" + row + ") that is not in the " +
                        "search window"
            )
        } else cellValues[colOffsets[col] + row - window.minJforI(col)] = value
    }

    override fun get(col: Int, row: Int): Double {
        return if (row < window.minJforI(col) || row > window.maxJforI(col)) OUT_OF_WINDOW_VALUE else cellValues[colOffsets[col] + row - window.minJforI(
            col
        )]
    }

    override fun size(): Int {
        return cellValues.size
    }

    companion object {
        // CONSTANTS
        private val OUT_OF_WINDOW_VALUE = Double.POSITIVE_INFINITY
    }

    // CONSTRUCTOR
    init {
        cellValues = DoubleArray(window.size())
        colOffsets = IntArray(window.maxI() + 1)

        // Fill in the offset matrix
        var currentOffset = 0
        for (i in window.minI()..window.maxI()) {
            colOffsets[i] = currentOffset
            currentOffset += window.maxJforI(i) - window.minJforI(i) + 1
        }
    } // end Constructor
} // end class MemoryResidentMatrix

class ParallelogramWindow(tsI: TimeSeries, tsJ: TimeSeries, searchRadius: Int) :
    SearchWindow(tsI.size(), tsJ.size()) {
    // CONSTRUCTOR
    init {

        // Find the coordinates of the parallelogram's corners..other than (minI,minJ) and (maxI, maxJ)
        val upperCornerI = Math.max(
            maxI() / 2.0 - searchRadius * (maxI().toDouble() / maxJ()),
            minI().toDouble()
        )
        val upperCornerJ = Math.min(
            maxJ() / 2.0 + searchRadius * (maxJ().toDouble() / maxI()),
            maxJ().toDouble()
        )
        val lowerCornerI = Math.min(
            maxI() / 2.0 + searchRadius * (maxI().toDouble() / maxJ()),
            maxI().toDouble()
        )
        val lowerCornerJ = Math.max(
            maxJ() / 2.0 - searchRadius * (maxJ().toDouble() / maxI()),
            minJ().toDouble()
        )

        // For each column determine the minimum and maximum row ranges that are in the paralellogram's window.
        for (i in 0 until tsI.size()) {
            val minJ: Int
            val maxJ: Int
            val isIlargest = tsI.size() >= tsJ.size()
            maxJ = if (i < upperCornerI) // left side of upper line
            {
                if (isIlargest) {
                    val interpRatio = i / upperCornerI
                    Math.round(interpRatio * upperCornerJ).toInt()
                } else {
                    val interpRatio = (i + 1) / upperCornerI
                    Math.round(interpRatio * upperCornerJ).toInt() - 1
                } // end if
            } else  // right side of upper line
            {
                if (isIlargest) {
                    val interpRatio = (i - upperCornerI) / (maxI() - upperCornerI)
                    Math.round(upperCornerJ + interpRatio * (maxJ() - upperCornerJ)).toInt()
                } else {
                    val interpRatio =
                        (i + 1 - upperCornerI) / (maxI() - upperCornerI)
                    Math.round(upperCornerJ + interpRatio * (maxJ() - upperCornerJ)).toInt() - 1
                } // end if
            } // end if
            minJ = if (i <= lowerCornerI) // left side of lower line
            {
                val interpRatio = i / lowerCornerI
                Math.round(interpRatio * lowerCornerJ).toInt()
            } else  // right side of lower line
            {
                val interpRatio = (i - lowerCornerI) / (maxI() - lowerCornerI)
                Math.round(lowerCornerJ + interpRatio * (maxJ() - lowerCornerJ)).toInt()
            } // end if
            super.markVisited(i, minJ)
            super.markVisited(i, maxJ)
        } // end for loop
    } // end Constructor
} // end class ParallelogramWindow

class TimeWarpInfo // CONSTRUCTOR
internal constructor(// PRIVATE DATA
    val distance: Double, val path: WarpPath
) {

    override fun toString(): String {
        return "(Warp Distance=$distance, Warp Path=$path)"
    }

} // end class TimeWarpInfo

internal class PartialWindowMatrix(private val window: SearchWindow) : CostMatrix {
    private var lastCol: DoubleArray
    private lateinit var currCol: DoubleArray
    private var currColIndex = 0
    private var minLastRow = 0
    private var minCurrRow: Int

    // PUBLIC FUNCTIONS
    override fun put(col: Int, row: Int, value: Double) {
        if (row < window.minJforI(col) || row > window.maxJforI(col)) {
            throw InternalError(
                "CostMatrix is filled in a cell (col=" + col + ", row=" + row + ") that is not in the " +
                        "search window"
            )
        } else {
            if (col == currColIndex) currCol[row - minCurrRow] =
                value else if (col == currColIndex - 1) lastCol[row - minLastRow] =
                value else if (col == currColIndex + 1) {
                lastCol = currCol
                minLastRow = minCurrRow
                currColIndex++
                currCol = DoubleArray(window.maxJforI(col) - window.minJforI(col) + 1)
                minCurrRow = window.minJforI(col)
                currCol[row - minCurrRow] = value
            } else throw InternalError("A PartialWindowMatrix can only fill in 2 adjacentcolumns at a time")
        } // end if
    } // end put(...)

    override fun get(col: Int, row: Int): Double {
        return if (row < window.minJforI(col) || row > window.maxJforI(col)) OUT_OF_WINDOW_VALUE else {
            if (col == currColIndex) currCol[row - minCurrRow] else if (col == currColIndex - 1) lastCol[row - minLastRow] else OUT_OF_WINDOW_VALUE
        } // end if
    } // end get(..)

    override fun size(): Int {
        return lastCol.size + currCol.size
    }

    fun windowSize(): Int {
        return window.size()
    }

    companion object {
        // PRIVATE DATA
        private val OUT_OF_WINDOW_VALUE = Double.POSITIVE_INFINITY
    }

    // CONSTRUCTOR
    init {
        if (window.maxI() > 0) {
            currCol = DoubleArray(window.maxJforI(1) - window.minJforI(1) + 1)
            currColIndex = 1
            minLastRow = window.minJforI(currColIndex - 1)
        } else currColIndex = 0
        minCurrRow = window.minJforI(currColIndex)
        lastCol = DoubleArray(window.maxJforI(0) - window.minJforI(0) + 1)
    } // end Constructor
} // end WindowMatrix


internal class SwapFileMatrix(// PRIVATE DATA
    private val window: SearchWindow
) : CostMatrix {

    // Private data needed to store the last 2 colums of the matrix.
    private var lastCol: DoubleArray
    private lateinit var currCol: DoubleArray
    private var currColIndex = 0
    private var minLastRow = 0
    private var minCurrRow: Int

    // Private data needed to read values from the swap file.
    private val swapFile: File
    private var cellValuesFile: RandomAccessFile? = null
    private val isSwapFileFreed: Boolean
    private val colOffsets: LongArray

    // PUBLIC FUNCTIONS
    override fun put(col: Int, row: Int, value: Double) {
        if (row < window.minJforI(col) || row > window.maxJforI(col)) {
            throw InternalError(
                "CostMatrix is filled in a cell (col=" + col + ", row=" + row + ") that is not in the " +
                        "search window"
            )
        } else {
            if (col == currColIndex) currCol[row - minCurrRow] =
                value else if (col == currColIndex - 1) {
                lastCol[row - minLastRow] = value
            } else if (col == currColIndex + 1) {
                // Write the last column to the swap file.
                try {
                    if (isSwapFileFreed) throw InternalError("The SwapFileMatrix has been freeded by the freeMem() method") else {
                        cellValuesFile!!.seek(cellValuesFile!!.length()) // move file poiter to end of file
                        colOffsets[currColIndex - 1] = cellValuesFile!!.filePointer

                        // Write an entire column to the swap file.
                        cellValuesFile!!.write(doubleArrayToByteArray(lastCol))
                    } // end if
                } catch (e: IOException) {
                    throw InternalError("Unable to fill the CostMatrix in the Swap file (IOException)")
                } // end try
                lastCol = currCol
                minLastRow = minCurrRow
                minCurrRow = window.minJforI(col)
                currColIndex++
                currCol = DoubleArray(window.maxJforI(col) - window.minJforI(col) + 1)
                currCol[row - minCurrRow] = value
            } else throw InternalError("A SwapFileMatrix can only fill in 2 adjacentcolumns at a time")
        } // end if
    } // end put(...)

    override fun get(col: Int, row: Int): Double {
        return if (row < window.minJforI(col) || row > window.maxJforI(col)) OUT_OF_WINDOW_VALUE else if (col == currColIndex) currCol[row - minCurrRow] else if (col == currColIndex - 1) lastCol[row - minLastRow] else {
            try {
                if (isSwapFileFreed) throw InternalError("The SwapFileMatrix has been freeded by the freeMem() method") else {
                    cellValuesFile!!.seek(colOffsets[col] + 8 * (row - window.minJforI(col)))
                    cellValuesFile!!.readDouble()
                } // end if
            } catch (e: IOException) {
                if (col > currColIndex) throw InternalError(
                    "The requested value is in the search window but has not been entered into " +
                            "the matrix: (col=" + col + "row=" + row + ")."
                ) else throw InternalError("Unable to read CostMatrix in the Swap file (IOException)")
            } // end try
        } // end if
    } // end get(..)

    // This method closes and delets the swap file when the object's finalize() mehtod is called.  This method will
    //    ONLY be called by the JVM if the object is garbage collected while the application is still running.
    //    This method must be called explicitly to guarantee that the swap file is deleted.
    @Throws(Throwable::class)


    protected fun finalize() {
        // Close and Delete the (possibly VERY large) swap file.
        try {
            if (!isSwapFileFreed) cellValuesFile!!.close()
        } catch (e: Exception) {
            System.err.println("unable to close swap file '" + swapFile.path + "' during finialization")
        } finally {
            swapFile.delete() // delete the swap file
            //super.finalize() // ensure that finalization continues for parent if an exception is thrown
        } // end try
    } // end finalize()

    override fun size(): Int {
        return window.size()
    }

    fun freeMem() {
        try {
            cellValuesFile!!.close()
        } catch (e: IOException) {
            System.err.println("unable to close swap file '" + swapFile.path + "'")
        } finally {
            if (!swapFile.delete()) System.err.println("unable to delete swap file '" + swapFile.path + "'")
        } // end try
    } // end freeMem

    companion object {
        // CONSTANTS
        private val OUT_OF_WINDOW_VALUE = Double.POSITIVE_INFINITY
        private val RAND_GEN = Random()
    }

    // CONSTRUCTOR
    init {
        if (window.maxI() > 0) {
            currCol = DoubleArray(window.maxJforI(1) - window.minJforI(1) + 1)
            currColIndex = 1
            minLastRow = window.minJforI(currColIndex - 1)
        } else  // special case for a <=1 point time series, less than 2 columns to fill in
            currColIndex = 0
        minCurrRow = window.minJforI(currColIndex)
        lastCol = DoubleArray(window.maxJforI(0) - window.minJforI(0) + 1)
        swapFile = File("swap" + RAND_GEN.nextLong())
        isSwapFileFreed = false
        //swapFile.deleteOnExit();
        colOffsets = LongArray(window.maxI() + 1)
        cellValuesFile = try {
            RandomAccessFile(swapFile, "rw")
        } catch (e: FileNotFoundException) {
            throw InternalError("ERROR:  Unable to create swap file: $swapFile")
        } // end try
    } // end Constructor
} // end class SwapFileMatrix


abstract class SearchWindow(tsIsize: Int, tsJsize: Int) {
    // PRIVATE DATA
    private val minValues: IntArray
    private val maxValues: IntArray
    private val maxJ: Int
    private var size: Int
    protected var modCount: Int
        private set

    // PUBLIC FUNCTIONS
    fun isInWindow(i: Int, j: Int): Boolean {
        return i >= minI() && i <= maxI() && minValues[i] <= j && maxValues[i] >= j
    }

    fun minI(): Int {
        return 0
    }

    fun maxI(): Int {
        return minValues.size - 1
    }

    fun minJ(): Int {
        return 0
    }

    fun maxJ(): Int {
        return maxJ
    }

    fun minJforI(i: Int): Int {
        return minValues[i]
    }

    fun maxJforI(i: Int): Int {
        return maxValues[i]
    }

    fun size(): Int {
        return size
    }

    // Iterates through all cells in the search window in the order that Dynamic
    //    Time Warping needs to evaluate them. (first to last column (0..maxI),
    //    bottom up  (o..maxJ))
    operator fun iterator(): Iterator<*> {
        return SearchWindowIterator(this)
    }

    override fun toString(): String {
        val outStr = StringBuffer()
        for (i in minI()..maxI()) {
            outStr.append("i=" + i + ", j=" + minValues[i] + "..." + maxValues[i])
            if (i != maxI()) outStr.append("\n")
        } // end for loop
        return outStr.toString()
    } // end toString()

    // PROTECTED FUNCTIONS
    //    Expands the current window by a s pecified radius.
    protected fun expandWindow(radius: Int) {
        if (radius > 0) {
            // Expand the search window by one before expanding by the remainder of the radius because the function
            //    "expandSearchWindow(.) may not work correctly if the path has a width of only 1.
            expandSearchWindow(1)
            expandSearchWindow(radius - 1)
        }
    }

    private fun expandSearchWindow(radius: Int) {
        if (radius > 0) // if radius <=0 then no search is necessary, use the current search window
        {
            // Add all cells in the current Window to an array, iterating through the window and expanding the window
            //    at the same time is not possible because the window can't be changed during iteration through the cells.
            val windowCells: ArrayList<Any?> = ArrayList<Any?>(size())
            val cellIter = this.iterator()
            while (cellIter.hasNext()) {
                windowCells.add(cellIter.next())
            }
            for (cell in windowCells.indices) {
                val currentCell = windowCells[cell] as ColMajorCell
                if (currentCell.col != minI() && currentCell.row != maxJ()) // move to upper left if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col - radius
                    val targetRow = currentCell.row + radius
                    if (targetCol >= minI() && targetRow <= maxJ()) markVisited(
                        targetCol,
                        targetRow
                    ) else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge =
                            Math.max(minI() - targetCol, targetRow - maxJ())
                        markVisited(targetCol + cellsPastEdge, targetRow - cellsPastEdge)
                    } // end if
                } // end if
                if (currentCell.row != maxJ()) // move up if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col
                    val targetRow = currentCell.row + radius
                    if (targetRow <= maxJ()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge = targetRow - maxJ()
                        markVisited(targetCol, targetRow - cellsPastEdge)
                    } // end if
                } // end if
                if (currentCell.col != maxI() && currentCell.row != maxJ()) // move to upper-right if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col + radius
                    val targetRow = currentCell.row + radius
                    if (targetCol <= maxI() && targetRow <= maxJ()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge =
                            Math.max(targetCol - maxI(), targetRow - maxJ())
                        markVisited(targetCol - cellsPastEdge, targetRow - cellsPastEdge)
                    } // end if
                } // end if
                if (currentCell.col != minI()) // move left if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col - radius
                    val targetRow = currentCell.row
                    if (targetCol >= minI()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge = minI() - targetCol
                        markVisited(targetCol + cellsPastEdge, targetRow)
                    } // end if
                } // end if
                if (currentCell.col != maxI()) // move right if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col + radius
                    val targetRow = currentCell.row
                    if (targetCol <= maxI()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge = targetCol - maxI()
                        markVisited(targetCol - cellsPastEdge, targetRow)
                    } // end if
                } // end if
                if (currentCell.col != minI() && currentCell.row != minJ()) // move to lower-left if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col - radius
                    val targetRow = currentCell.row - radius
                    if (targetCol >= minI() && targetRow >= minJ()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge =
                            Math.max(minI() - targetCol, minJ() - targetRow)
                        markVisited(targetCol + cellsPastEdge, targetRow + cellsPastEdge)
                    } // end if
                } // end if
                if (currentCell.row != minJ()) // move down if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col
                    val targetRow = currentCell.row - radius
                    if (targetRow >= minJ()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge = minJ() - targetRow
                        markVisited(targetCol, targetRow + cellsPastEdge)
                    } // end if
                } // end if
                if (currentCell.col != maxI() && currentCell.row != minJ()) // move to lower-right if possible
                {
                    // Either extend full search radius or some fraction until edges of matrix are met.
                    val targetCol = currentCell.col + radius
                    val targetRow = currentCell.row - radius
                    if (targetCol <= maxI() && targetRow >= minJ()) markVisited(
                        targetCol,
                        targetRow
                    ) // radius does not go past the edges of the matrix
                    else {
                        // Expand the window only to the edge of the matrix.
                        val cellsPastEdge =
                            Math.max(targetCol - maxI(), minJ() - targetRow)
                        markVisited(targetCol - cellsPastEdge, targetRow + cellsPastEdge)
                    } // end if
                } // end if
            } // end for loop
        } // end if
    } // end expandWindow(.)

    // Raturns true if the window is modified.
    protected fun markVisited(col: Int, row: Int) {
        if (minValues[col] == -1) // first value is entered in the column
        {
            minValues[col] = row
            maxValues[col] = row
            size++
            modCount++ // stucture has been changed
            //return true;
        } else if (minValues[col] > row) // minimum range in the column is expanded
        {
            size += minValues[col] - row
            minValues[col] = row
            modCount++ // stucture has been changed
        } else if (maxValues[col] < row) // maximum range in the column is expanded
        {
            size += row - maxValues[col]
            maxValues[col] = row
            modCount++
        } // end if
    } // end markVisited(.)

    // A private class that is a fail-fast iterator through the search window.
    inner class SearchWindowIterator  constructor( val window: SearchWindow) :
        MutableIterator<Any?> {
        // PRIVATE DATA
        private var currentI: Int
        private var currentJ: Int
        private var hasMoreElements: Boolean
        private val expectedModCount: Int

        // PUBLIC FUNCTIONS
        override fun hasNext(): Boolean {
            return hasMoreElements
        }

        override fun next(): Any {
            return if (modCount != expectedModCount) throw ConcurrentModificationException() else if (!hasMoreElements) throw NoSuchElementException() else {
                val cell = ColMajorCell(currentI, currentJ)
                if (++currentJ > window.maxJforI(currentI)) {
                    if (++currentI <= window.maxI()) currentJ =
                        window.minJforI(currentI) else hasMoreElements = false
                } // end if
                cell
            } // end if
        } // end next()

        override fun remove() {
            throw UnsupportedOperationException()
        }

        // CONSTRUCTOR
        init {
            // Intiialize values
            hasMoreElements = window.size() > 0
            currentI = window.minI()
            currentJ = window.minJ()
            expectedModCount = window.modCount
        }
    } // end inner class SearchWindowIterator

    // CONSTRUCTOR
    init {
        minValues = IntArray(tsIsize)
        maxValues = IntArray(tsIsize)
        Arrays.fill(minValues, -1)
        maxJ = tsJsize - 1
        size = 0
        modCount = 0
    }
} // end SearchWindow


class Array2D {
    // PRIVATE DATA
    private val rows // ArrayList of ArrayList (an array of rows in the array)
            : ArrayList<Any?>
    private var numOfElements: Int

    // CONSTRUCTOR
    constructor() {
        rows = ArrayList<Any?>()
        numOfElements = 0
    }

    constructor(array: Array2D) {
        rows = ArrayList(array.rows)
        numOfElements = array.numOfElements
    }

    // PUBLIC FU?NCTIONS
    fun clear() {
        rows.clear()
        numOfElements = 0
    }

    fun size(): Int {
        return numOfElements
    }

    fun numOfRows(): Int {
        return rows.size
    }

    fun getSizeOfRow(row: Int): Int {
        return (rows[row] as ArrayList<Any?>).size
    }

    operator fun get(row: Int, col: Int): Any {
        return (rows[row] as ArrayList<*>)[col]
    }

    operator fun set(row: Int, col: Int, newVal: Any?) {
        (rows[row] as ArrayList<Any?>).set(col, newVal)
    }

    fun addToEndOfRow(row: Int, value: Any?) {
        (rows[row] as ArrayList<Any?>).add(value)
        numOfElements++
    }

    fun addAllToEndOfRow(row: Int, objects: Collection<*>) {
        val i = objects.iterator()
        while (i.hasNext()) {
            (rows[row] as ArrayList<Any?>).add(i.next())
            numOfElements++
        }
    }

    fun addToNewFirstRow(value: Any?) {
        val newFirstRow: ArrayList<Any?> = ArrayList<Any?>(1)
        newFirstRow.add(value)
        rows.add(0, newFirstRow)
        numOfElements++
    }

    fun addToNewLastRow(value: Any?) {
        val newLastRow: ArrayList<Any?> = ArrayList<Any?>(1)
        newLastRow.add(value)
        rows.add(newLastRow)
        numOfElements++
    }

    fun addAllToNewLastRow(objects: Collection<*>) {
        val i = objects.iterator()
        val newLastRow: ArrayList<Any?> = ArrayList<Any?>(1)
        while (i.hasNext()) {
            newLastRow.add(i.next())
            numOfElements++
        }
        rows.add(newLastRow)
    }

    fun removeFirstRow() {
        numOfElements -= (rows[0] as ArrayList<Any?>).size
        rows.removeAt(0)
    }

    fun removeLastRow() {
        numOfElements -= (rows[rows.size - 1] as ArrayList<Any?>).size
        rows.removeAt(rows.size - 1)
    }

    override fun toString(): String {
        var outStr = ""
        for (r in rows.indices) {
            val currentRow = rows[r] as ArrayList<*>
            for (c in currentRow.indices) {
                outStr += currentRow[c]
                outStr += if (c == currentRow.size - 1) "\n" else ","
            }
        } // end for
        return outStr
    }
} // end class matrix.Array2D

internal interface CostMatrix {
    fun put(col: Int, row: Int, value: Double)
    operator fun get(col: Int, row: Int): Double
    fun size(): Int
} // end interface CostMatrix


class ColMajorCell // end Constructor
    (val col: Int, val row: Int) {

    override fun equals(o: Any?): Boolean {
        return o is ColMajorCell &&
                o.col == col &&
                o.row == row
    }

    override fun hashCode(): Int {
        return (1 shl col) + row
    }

    override fun toString(): String {
        return "($col,$row)"
    }

    fun toInt(): Int {
        return col.toInt()
    }

} // end class ColMajorCell


class TimeSeriesPoint {
    // PRIVATE DATA
    private var measurements: DoubleArray
    private var hashCode: Int

    // CONSTRUCTORS
    constructor(values: DoubleArray) {
        hashCode = 0
        measurements = DoubleArray(values.size)
        for (x in values.indices) {
            hashCode += values[x].hashCode()
            measurements[x] = values[x]
        }
    }

    constructor(values: Collection<*>) {
        measurements = DoubleArray(values.size)
        hashCode = 0
        val i = values.iterator()
        var index = 0
        while (i.hasNext()) {
            val nextElement = i.next()!!
            if (nextElement is Double) measurements[index] =
                nextElement.toDouble() else if (nextElement is Int) measurements[index] =
                nextElement.toDouble() else if (nextElement is BigInteger) measurements[index] =
                nextElement.toDouble() else throw InternalError(
                "ERROR:  The element " + nextElement +
                        " is not a valid numeric type"
            )
            hashCode += measurements[index].hashCode()
            index++
        } // end while loop
    } // end constructor

    // FUNCTIONS
    operator fun get(dimension: Int): Double {
        return measurements[dimension]
    }

    operator fun set(dimension: Int, newValue: Double) {
        hashCode -= measurements[dimension].hashCode()
        measurements[dimension] = newValue
        hashCode += newValue.hashCode()
    }

    fun toArray(): DoubleArray {
        return measurements
    }

    fun size(): Int {
        return measurements.size
    }

    override fun toString(): String {
        var outStr = "("
        for (x in measurements.indices) {
            outStr += measurements[x]
            if (x < measurements.size - 1) outStr += ","
        }
        outStr += ")"
        return outStr
    } // end toString()

    override fun equals(o: Any?): Boolean {
        return if (this === o) true else if (o is TimeSeriesPoint) {
            val testValues = o.toArray()
            if (testValues.size == measurements.size) {
                for (x in measurements.indices) if (measurements[x] != testValues[x]) return false
                true
            } else false
        } else false
    } // end public boolean equals

    override fun hashCode(): Int {
        return hashCode
    }
} // end class TimeSeriesPoint


object TypeConversions {
    fun doubleToByteArray(number: Double): ByteArray {
        // double to long representation
        val longNum = java.lang.Double.doubleToLongBits(number)

        // long to 8 bytes
        return byteArrayOf(
            (longNum ushr 56 and 0xFF).toByte(),
            (longNum ushr 48 and 0xFF).toByte(),
            (longNum ushr 40 and 0xFF).toByte(),
            (longNum ushr 32 and 0xFF).toByte(),
            (longNum ushr 24 and 0xFF).toByte(),
            (longNum ushr 16 and 0xFF).toByte(),
            (longNum ushr 8 and 0xFF).toByte(),
            (longNum ushr 0 and 0xFF).toByte()
        )
    } // end doubleToByte(.)

    fun doubleArrayToByteArray(numbers: DoubleArray): ByteArray {
        val doubleSize = 8 // 8 byes in a double
        val byteArray = ByteArray(numbers.size * doubleSize)
        for (x in numbers.indices) System.arraycopy(
            doubleToByteArray(
                numbers[x]
            ), 0, byteArray, x * doubleSize, doubleSize
        )
        return byteArray
    } // end doubleArrayToByteArray(.)
} // end class Typeconversions


/**
 * This class...
 *
 * @author Stan Salvador, stansalvador@hotmail.com
 * @version last changed: Jun 30, 2004
 * @see
 * @since Jun 30, 2004
 */
class SineWave(length: Int, cycles: Double, noise: Double) :
    TimeSeries(1) {
    companion object {
        private val rand = Random()
    }

    // CONSTRUCTORS
    init {

        // final Random rand = new Random();
        for (x in 0 until length) {
            val nextPoint =
                Math.sin(x.toDouble() / length * 2.0 * Math.PI * cycles) + rand.nextGaussian() * noise
            super.addLast(x.toDouble(), TimeSeriesPoint(doubleArrayOf(nextPoint)))
        }
    } // PUBLIC FUNCTIONS
} // end class SineWave


open class TimeSeries {
    // PRIVATE DATA
    private val labels // labels for each column
            : ArrayList<Any?>
    private val timeReadings // ArrayList of Double
            : ArrayList<Any?>
    private val tsArray // ArrayList of TimeSeriesPoint.. no time
            : ArrayList<Any?>

    /**
     * Ad-hoc constructor from a double array
     * ONLY ONE FIXED COLUMN
     */
    constructor(x: DoubleArray) {
        // DS to create labels and values
        val counter =
            AtomicInteger()
        var tsValues: TimeSeriesPoint?

        // Only one time series
        timeReadings = ArrayList<Any?>()
        var values: ArrayList<Any?>

        // Each value is one incremented label
        labels = ArrayList<Any?>()
        labels.add("Time")
        labels.add("c1")
        tsArray = ArrayList<Any?>()
        for (value in x) {
            timeReadings.add(counter.getAndIncrement().toDouble())
            values = ArrayList<Any?>()
            values.add(value)
            tsValues = TimeSeriesPoint(values)
            tsArray.add(tsValues)
        }
    }

    // TODO don't use defaults delimiter/1stColTime... determine if not specified
    // CONSTRUCTORS                                                                   // TODO method to peek at determined delimiter, 1st col time
    internal constructor() {
        labels = ArrayList<Any?>() // TODO isLabeled constuctor options?
        timeReadings = ArrayList<Any?>()
        tsArray = ArrayList<Any?>()
    }

    constructor(numOfDimensions: Int) : this() {
        labels.add("Time")
        for (x in 0 until numOfDimensions) labels.add("" + x)
    }

    // Copy Constructor
    constructor(origTS: TimeSeries) {
        labels = ArrayList(origTS.labels)
        timeReadings = ArrayList(origTS.timeReadings)
        tsArray = ArrayList(origTS.tsArray)
    }

    constructor(inputFile: String, isFirstColTime: Boolean) : this(
        inputFile,
        ZERO_ARRAY,
        isFirstColTime
    ) {
    }

    constructor(inputFile: String, delimiter: Char) : this(
        inputFile,
        ZERO_ARRAY,
        DEFAULT_IS_TIME_1ST_COL,
        DEFAULT_IS_LABELED,
        delimiter
    ) {
    }

    constructor(
        inputFile: String,
        isFirstColTime: Boolean,
        delimiter: Char
    ) : this(
        inputFile,
        ZERO_ARRAY,
        isFirstColTime,
        DEFAULT_IS_LABELED,
        delimiter
    ) {
    }

    constructor(
        inputFile: String,
        isFirstColTime: Boolean,
        isLabeled: Boolean,
        delimiter: Char
    ) : this(inputFile, ZERO_ARRAY, isFirstColTime, isLabeled, delimiter) {
    }

    @JvmOverloads
    constructor(
        inputFile: String,
        colToInclude: IntArray?,
        isFirstColTime: Boolean,
        isLabeled: Boolean = DEFAULT_IS_LABELED,
        delimiter: Char = DEFAULT_DELIMITER
    ) : this() {
        try {
            // Record the Label names (fropm the top row.of the input file).
            var br =
                BufferedReader(FileReader(inputFile)) // open the input file
            var line = br.readLine() // the top row that contains attribiute names.
            var st =
                StringTokenizer(line, delimiter.toString())
            if (isLabeled) {
                var currentCol = 0
                while (st.hasMoreTokens()) {
                    val currentToken = st.nextToken()
                    if (colToInclude!!.size == 0 || contains(
                            colToInclude,
                            currentCol
                        )
                    ) labels.add(currentToken)
                    currentCol++
                } // end while loop

                // Make sure that the first column is labeled is for Time.
                if (labels.size == 0) throw InternalError(
                    "ERROR:  The first row must contain label " +
                            "information, it is empty!"
                ) else if (!isFirstColTime) labels.add(
                    0,
                    "Time"
                ) else if (isFirstColTime && !(labels[0] as String).equals(
                        "Time",
                        ignoreCase = true
                    )
                ) throw InternalError(
                    "ERROR:  The time column (1st col) in a time series must be labeled as 'Time', '" +
                            labels[0] + "' was found instead"
                )
            } else  // time series file is not labeled
            {
                if (colToInclude == null || colToInclude.size == 0) {
                    labels.add("Time")
                    if (isFirstColTime) st.nextToken()
                    var currentCol = 1 // TODO input fails gracefully
                    while (st.hasMoreTokens()) {
                        st.nextToken()
                        labels.add(("c" + currentCol++)) // TODO add measurement with no time
                    }
                } else {
                    Arrays.sort(colToInclude)
                    labels.add("Time")
                    for (c in colToInclude.indices) if (colToInclude[c] > 0) labels.add(("c$c")) // TODO change to letterNum
                } // end if

                // Close and re-open the file.
                br.close()
                br = BufferedReader(FileReader(inputFile)) // open the input file
            } // end if


            // Read in all of the values in the data file.
            while (br.readLine().also { line = it } != null) // read lines until end of file
            {
                if (line.length > 0) // ignore empty lines
                {
                    st = StringTokenizer(line, delimiter.toString())

                    // Make sure that the current line has the correct number of
                    //    currentLineValues in it.
                    //           if (st.countTokens() != (labels.size()+ignoredCol))
                    //              throw new InternalError("ERROR:  Line " + (tsArray.size()+1) +
                    //                                      "contains the wrong number of currentLineValues. " +
                    //                                      "expected:  " + (labels.size()+ignoredCol) + ", " +
                    //                                      "found: " + st.countTokens());

                    // Read all currentLineValues in the current line
                    val currentLineValues: ArrayList<Any?> =
                        ArrayList<Any?>()
                    var currentCol = 0
                    while (st.hasMoreTokens()) {
                        val currentToken = st.nextToken()
                        if (colToInclude!!.size == 0 || contains(
                                colToInclude,
                                currentCol
                            )
                        ) {
                            // Attempt to parse the next value to a Double value.
                            val nextValue: Double
                            nextValue = try {
                                java.lang.Double.valueOf(currentToken)
                            } catch (e: NumberFormatException) {
                                throw InternalError("ERROR:  '$currentToken' is not a valid number")
                            }
                            currentLineValues.add(nextValue)
                        } // end if
                        currentCol++
                    } // end while loop

                    // Update the private data with the current Row that has been
                    //    read.
                    if (isFirstColTime) timeReadings.add(currentLineValues[0]) else timeReadings.add(
                        timeReadings.size
                    )
                    val firstMeasurement: Int
                    firstMeasurement = if (isFirstColTime) 1 else 0
                    val readings = TimeSeriesPoint(
                        currentLineValues.subList(
                            firstMeasurement,
                            currentLineValues.size
                        )
                    )
                    tsArray.add(readings)
                    //               timeValueMap.put(timeReadings.get(timeReadings.size()-1), readings);
                } // end if
            } // end while loop
            br.close()
        } catch (e: FileNotFoundException) {
            throw InternalError("ERROR:  The file '$inputFile' was not found.")
        } catch (e: IOException) {
            throw InternalError("ERROR:  Problem reading the file '$inputFile'.")
        } // end try block
    } // end constructor

    // FUNCTIONS
    @Throws(IOException::class)
    fun save(outFile: File?) {
        val out = PrintWriter(FileOutputStream(outFile))
        out.write(this.toString())
        out.flush()
        out.close()
    }

    fun clear() {
        labels.clear()
        timeReadings.clear()
        tsArray.clear()
        //     timeValueMap.clear();
    }

    fun size(): Int {
        return timeReadings.size
    }

    fun numOfPts(): Int {
        return size()
    }

    fun numOfDimensions(): Int {
        return labels.size - 1
    }

    fun getTimeAtNthPoint(n: Int): Double {
        return (timeReadings[n] as Double).toDouble()
    }

    fun getLabel(index: Int): String {
        return labels[index] as String
    }

    val labelsArr: Array<String?>
        get() {
            val labelArr = arrayOfNulls<String>(labels.size)
            for (x in labels.indices) labelArr[x] = labels[x] as String
            return labelArr
        }

    fun getLabels(): ArrayList<Any?> {
        return labels
    }

    fun setLabels(newLabels: ArrayList<Any?>) {
        labels.clear()
        for (x in newLabels.indices) labels.add(newLabels[x])
    }

    fun setLabels(newLabels: Array<String?>) {
        labels.clear()
        for (x in newLabels.indices) labels.add(newLabels[x])
    }

    fun getMeasurement(pointIndex: Int, valueIndex: Int): Double {
        return (tsArray[pointIndex] as TimeSeriesPoint).get(valueIndex)
    }

    fun getMeasurement(pointIndex: Int, valueLabel: String): Double {
        val valueIndex = labels.indexOf(valueLabel)
        if (valueIndex < 0) throw InternalError(
            "ERROR:  the label '" + valueLabel + "' was " +
                    "not one of:  " + labels
        )
        return (tsArray[pointIndex] as TimeSeriesPoint).get(valueIndex - 1)
    }

    fun getMeasurementVector(pointIndex: Int): DoubleArray {
        return (tsArray[pointIndex] as TimeSeriesPoint).toArray()
    }

    fun setMeasurement(pointIndex: Int, valueIndex: Int, newValue: Double) {
        (tsArray[pointIndex] as TimeSeriesPoint).set(valueIndex, newValue)
    }

    fun addFirst(time: Double, values: TimeSeriesPoint) {
        if (labels.size != values.size() + 1) throw InternalError(
            "ERROR:  The TimeSeriesPoint: " + values +
                    " contains the wrong number of values. " +
                    "expected:  " + labels.size + ", " +
                    "found: " + values.size()
        )
        if (time >= (timeReadings[0] as Double).toDouble()) throw InternalError(
            "ERROR:  The point being inserted into the " +
                    "beginning of the time series does not have " +
                    "the correct time sequence. "
        )
        timeReadings.add(0, time)
        tsArray.add(0, values)
    } // end addFirst(..)

    fun addLast(time: Double, values: TimeSeriesPoint) {
        if (labels.size != values.size() + 1) throw InternalError(
            "ERROR:  The TimeSeriesPoint: " + values +
                    " contains the wrong number of values. " +
                    "expected:  " + labels.size + ", " +
                    "found: " + values.size()
        )
        if (size() > 0 && time <= (timeReadings[timeReadings.size - 1] as Double).toDouble()) throw InternalError(
            "ERROR:  The point being inserted at the " +
                    "end of the time series does not have " +
                    "the correct time sequence. "
        )
        timeReadings.add(time)
        tsArray.add(values)
    } // end addLast(..)

    fun removeFirst() {
        if (size() == 0) System.err.println("WARNING:  TimeSeriesPoint:removeFirst() called on an empty time series!") else {
            timeReadings.removeAt(0)
            tsArray.removeAt(0)
        } // end if
    } // end removeFirst()

    fun removeLast() {
        if (size() == 0) System.err.println("WARNING:  TimeSeriesPoint:removeLast() called on an empty time series!") else {
            tsArray.removeAt(timeReadings.size - 1)
            timeReadings.removeAt(timeReadings.size - 1)
        } // end if
    } // end removeFirst()

    fun normalize() {
        // Calculate the mean of each FD.
        val mean = DoubleArray(numOfDimensions())
        for (col in 0 until numOfDimensions()) {
            var currentSum = 0.0
            for (row in 0 until size()) currentSum += this.getMeasurement(row, col)
            mean[col] = currentSum / size()
        } // end for loop

        // Calculate the standard deviation of each FD.
        val stdDev = DoubleArray(numOfDimensions())
        for (col in 0 until numOfDimensions()) {
            var variance = 0.0
            for (row in 0 until size()) variance += Math.abs(
                getMeasurement(
                    row,
                    col
                ) - mean[col]
            )
            stdDev[col] = variance / size()
        } // end for loop


        // Normalize the values in the data using the mean and standard deviation
        //    for each FD.  =>  Xrc = (Xrc-Mc)/SDc
        for (row in 0 until size()) {
            for (col in 0 until numOfDimensions()) {
                // Normalize data point.
                if (stdDev[col] == 0.0) // prevent divide by zero errors
                    setMeasurement(row, col, 0.0) // stdDev is zero means all pts identical
                else  // typical case
                    setMeasurement(
                        row,
                        col,
                        (getMeasurement(row, col) - mean[col]) / stdDev[col]
                    )
            } // end for loop
        } // end for loop
    } // end normalize();

    override fun toString(): String {
        val outStr = StringBuffer()
        /*
      // Write labels
      for (int x=0; x<labels.size(); x++)
      {
         outStr.append(labels.get(x));
         if (x < labels.size()-1)
            outStr.append(",");
         else
            outStr.append("\n");
      }  // end for loop
*/
        // Write the data for each row.
        for (r in timeReadings.indices) {
            // Time
//         outStr.append(timeReadings.get(r).toString());

            // The rest of the value on the row.
            val values: TimeSeriesPoint = tsArray[r] as TimeSeriesPoint
            for (c in 0 until values.size()) outStr.append(values.get(c))
            if (r < timeReadings.size - 1) outStr.append("\n")
        } // end for loop
        return outStr.toString()
    } // end toString()

    protected fun setMaxCapacity(capacity: Int) {
        timeReadings.ensureCapacity(capacity)
        tsArray.ensureCapacity(capacity)
    }

    companion object {
        private val ZERO_ARRAY = IntArray(0)
        private const val DEFAULT_IS_TIME_1ST_COL = true
        private const val DEFAULT_DELIMITER = ','
        private const val DEFAULT_IS_LABELED = true

        // Returns the first non-digit (and not a '.') character in a file under the
        //    assumption that it is the delimiter in the file.
        private fun determineDelimiter(filePath: String): Char {
            val DEFAULT_DELIMITER = ','
            return try {
                val `in` =
                    BufferedReader(FileReader(filePath))
                var line = `in`.readLine().trim { it <= ' ' } // read first line
                if (!Character.isDigit(line[0])) // go to 2nd line if 1st line appears to be labels
                    line = `in`.readLine()
                `in`.close()

                // Searches the 2nd line of the file until a non-number character is
                //    found.  The delimiter is assumed to be that character.
                //    numbers, minus signs, periods, and 'E' (exponent) are accepted
                //    number characters.
                for (x in 0 until line.length) {
                    if (!Character.isDigit(line[x]) && line[x] != '.' && line[x] != '-' &&
                        Character.toUpperCase(line[x]) != 'E'
                    ) return line[x]
                }

                // No delimiters were found, which means that there must be only one column
                //    A delimiter does not need to be known to read this file.
                DEFAULT_DELIMITER
            } catch (e: IOException) {
                DEFAULT_DELIMITER
            }
        } // end determineDelimiter(.)

        private fun extractFirstNumber(str: String): Double {
            val numStr = StringBuffer()

            // Keep adding characters onto numStr until a non-number character
            //    is reached.
            for (x in 0 until str.length) {
                if (Character.isDigit(str[x]) || str[x] == '.' || str[x] == '-' ||
                    Character.toUpperCase(str[x]) == 'E'
                ) numStr.append(str[x]) else numStr.toString().toDouble()
            } // end for loop
            return (-1).toDouble()
        }

        // Automatically determines if the first column in a file is time measurements.
        //    It assumes that a column of time will have equal spacing between all
        //    values.
        private fun determineIsFirstColTime(filePath: String): Boolean {
            val DEFAULT_VALUE = false
            return try {
                val `in` =
                    BufferedReader(FileReader(filePath))

                // This parameter is the percentage of flexibility that is permitted from
                //    a perfectly even distribution of time values.  This function will
                //    most likely not work if this is set to zero because of floating-
                //    point math roundoff errors.
                //    (a setting of '0.05' is '5 percent')
                val EQUALITY_FLEXIBILITY_PCT = 0.001
                val NUM_OF_VALUES_TO_CMP = 100 // $ of time values to look examine
                val possibleTimeValues: Vector<Any?> =
                    Vector<Any?>(NUM_OF_VALUES_TO_CMP) // 'stores numOfValuesToCompare' values

                // Read the first 'numOfValuesToCompare' possible time values from the file
                //    and store them in 'possibleTimeValues'.
                var line = `in`.readLine()
                while (possibleTimeValues.size < NUM_OF_VALUES_TO_CMP && `in`.readLine()
                        .also { line = it } != null
                ) possibleTimeValues.add(extractFirstNumber(line))
                if (possibleTimeValues.size <= 1) return DEFAULT_VALUE

                // See if there is equal spacing (with a flexibility of
                //    'equalityFlexibilityFactor') between all values in              // TODO TimeSeries is now messy...in need of design
                //    'possibleTimeValues'.
                if (possibleTimeValues.size > 1 && possibleTimeValues[1] == possibleTimeValues[0]
                ) return DEFAULT_VALUE // special case needed for very flat data
                val expectedDiff =
                    (possibleTimeValues[1] as Double).toDouble() -
                            (possibleTimeValues[0] as Double).toDouble()
                val flexibility = expectedDiff * EQUALITY_FLEXIBILITY_PCT
                for (x in 1 until possibleTimeValues.size) {
                    if (Math.abs(
                            (possibleTimeValues[x] as Double).toDouble() -
                                    (possibleTimeValues[x - 1] as Double).toDouble() - expectedDiff
                        )
                        > Math.abs(flexibility)
                    ) {
                        return false
                    }
                } // end for loop
                true
            } catch (e: IOException) {
                DEFAULT_VALUE
            }
        } // end determineIsFirstColTime(.)
    }
} // end class TimeSeries


/**
 * This class...
 *
 * @author Stan Salvador, stansalvador@hotmail.com
 * @version last changed: 06/01/2004
 * @see
 * @since 06/01/2004
 */
object Arrays {
    fun toPrimitiveArray(objArr: Array<Int>): IntArray {
        val primArr = IntArray(objArr.size)
        for (x in objArr.indices) primArr[x] = objArr[x]
        return primArr
    }

    fun toIntArray(c: Collection<*>): IntArray {
        return toPrimitiveArray(c.toTypedArray() as Array<Int>)
    }

    fun toCollection(arr: BooleanArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: ByteArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: CharArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: DoubleArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: FloatArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: IntArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: LongArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: ShortArray): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add(arr[x])
        return collection
    }

    fun toCollection(arr: Array<String?>): Collection<*> {
        val collection: MutableList<Any?> = ArrayList<Any?>(arr.size)
        for (x in arr.indices) collection.add((arr[x]))
        return collection
    }

    fun contains(arr: BooleanArray, `val`: Boolean): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: ByteArray, `val`: Byte): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: CharArray, `val`: Char): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: DoubleArray, `val`: Double): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: FloatArray, `val`: Float): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: IntArray, `val`: Int): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: LongArray, `val`: Long): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: ShortArray, `val`: Short): Boolean {
        for (x in arr.indices) if (arr[x] == `val`) return true
        return false
    }

    fun contains(arr: Array<String>, `val`: String): Boolean {
        for (x in arr.indices) if (arr[x] === `val`) return true
        return false
    }
} // end class Arrays




interface DistanceFunction {
    fun calcDistance(vector1: DoubleArray?, vector2: DoubleArray?): Double
}

object DistanceFunctionFactory {
    var EUCLIDEAN_DIST_FN: DistanceFunction = EuclideanDistance()
    var MANHATTAN_DIST_FN: DistanceFunction = ManhattanDistance()
    var BINARY_DIST_FN: DistanceFunction = BinaryDistance()
    fun getDistFnByName(distFnName: String): DistanceFunction {
        return if (distFnName == "EuclideanDistance") {
            EUCLIDEAN_DIST_FN
        } else if (distFnName == "ManhattanDistance") {
            MANHATTAN_DIST_FN
        } else if (distFnName == "BinaryDistance") {
            BINARY_DIST_FN
        } else {
            throw IllegalArgumentException("There is no DistanceFunction for the name $distFnName")
        } // end if
    }
}

class BinaryDistance : DistanceFunction {
    override fun calcDistance(
        vector1: DoubleArray?,
        vector2: DoubleArray?
    ): Double {
        if (vector1!!.size != vector2!!.size) throw InternalError(
            "ERROR:  cannot calculate the distance "
                    + "between vectors of different sizes."
        ) else if (Arrays.equals(vector1, vector2)) {
            return 0.0
        } else {
            return 1.0
        }
    }
}


class EuclideanDistance : DistanceFunction {
    override fun calcDistance(
        vector1: DoubleArray?,
        vector2: DoubleArray?
    ): Double {
        if (vector1!!.size != vector2!!.size) throw InternalError(
            "ERROR:  cannot calculate the distance "
                    + "between vectors of different sizes."
        )
        var sqSum = 0.0
        for (x in vector1.indices) sqSum += Math.pow(vector1[x] - vector2[x], 2.0)
        return Math.sqrt(sqSum)
    } // end class euclideanDist(..)
}

class ManhattanDistance : DistanceFunction {
    override fun calcDistance(
        vector1: DoubleArray?,
        vector2: DoubleArray?
    ): Double {
        if (vector1!!.size != vector2!!.size) throw InternalError(
            "ERROR:  cannot calculate the distance "
                    + "between vectors of different sizes."
        )
        var diffSum = 0.0
        for (x in vector1.indices) diffSum += Math.abs(vector1[x] - vector2[x])
        return diffSum
    }
}




class ExpandedResWindow(
    tsI: TimeSeries, tsJ: TimeSeries, shrunkI: PAA, shrunkJ: PAA,
    shrunkWarpPath: WarpPath, searchRadius: Int
) :
    SearchWindow(tsI.size(), tsJ.size()) {
    init {
        // Initialize the private data in the super class.

        // Variables to keep track of the current location of the higher resolution projected path.
        var currentI: Int = shrunkWarpPath.minI()
        var currentJ: Int = shrunkWarpPath.minJ()

        // Variables to keep track of the last part of the low-resolution warp path that was evaluated
        //    (to determine direction).
        var lastWarpedI = Int.MAX_VALUE
        var lastWarpedJ = Int.MAX_VALUE

        // For each part of the low-resolution warp path, project that path to the higher resolution by filling in the
        //    path's corresponding cells at the higher resolution.
        for (w in 0 until shrunkWarpPath.size()) {
            val currentCell: ColMajorCell = shrunkWarpPath.get(w)
            val warpedI = currentCell.col
            val warpedJ = currentCell.row
            val blockISize = shrunkI.aggregatePtSize(warpedI)
            val blockJSize = shrunkJ.aggregatePtSize(warpedJ)

            // If the path moved up or diagonally, then the next cell's values on the J axis will be larger.
            if (warpedJ > lastWarpedJ) currentJ += shrunkJ.aggregatePtSize(lastWarpedJ)

            // If the path moved up or diagonally, then the next cell's values on the J axis will be larger.
            if (warpedI > lastWarpedI) currentI += shrunkI.aggregatePtSize(lastWarpedI)

            // If a diagonal move was performed, add 2 cells to the edges of the 2 blocks in the projected path to create
            //    a continuous path (path with even width...avoid a path of boxes connected only at their corners).
            //                        |_|_|x|x|     then mark      |_|_|x|x|
            //    ex: projected path: |_|_|x|x|  --2 more cells->  |_|X|x|x|
            //                        |x|x|_|_|        (X's)       |x|x|X|_|
            //                        |x|x|_|_|                    |x|x|_|_|
            if (warpedJ > lastWarpedJ && warpedI > lastWarpedI) {
                super.markVisited(currentI - 1, currentJ)
                super.markVisited(currentI, currentJ - 1)
            } // end if

            // Fill in the cells that are created by a projection from the cell in the low-resolution warp path to a
            //    higher resolution.
            for (x in 0 until blockISize) {
                super.markVisited(currentI + x, currentJ)
                super.markVisited(currentI + x, currentJ + blockJSize - 1)
            } // end for loop

            // Record the last position in the warp path so the direction of the path can be determined when the next
            //    position of the path is evaluated.
            lastWarpedI = warpedI
            lastWarpedJ = warpedJ
        } // end for loop


        // Expand the size of the projected warp path by the specified width.
        super.expandWindow(searchRadius)
    } // end Constructor
} // end class ExpandedResWindow

internal class WindowMatrix(searchWindow: SearchWindow) : CostMatrix {
    // PRIVATE DATA
    private var windowCells: CostMatrix? = null

    // PUBLIC FUNCTIONS
    override fun put(col: Int, row: Int, value: Double) {
        windowCells!!.put(col, row, value)
    }

    override fun get(col: Int, row: Int): Double {
        return windowCells!![col, row]
    }

    override fun size(): Int {
        return windowCells!!.size()
    }

    fun freeMem() {
        // Resources only need to be freed for a SwapFileMatrix.
        if (windowCells is SwapFileMatrix) {
            try {
                (windowCells as SwapFileMatrix).freeMem()
            } catch (t: Throwable) {
                // ignore the exception
            } // end try
        } // end if
    } // end freeMem()

    // CONSTRUCTOR
    init {
        try {
            windowCells = MemoryResidentMatrix(searchWindow)
        } catch (e: OutOfMemoryError) {
            System.err.println(
                "Ran out of memory initializing window matrix, all cells in the window cannot fit into " +
                        "main memory.  Will use a swap file instead (will run ~50% slower)"
            )
            System.gc()
            windowCells = SwapFileMatrix(searchWindow)
        } // end try
    } // end Constructor
} // end WindowMatrix

/**
 * This class...
 *
 * @author Stan Salvador, stansalvador@hotmail.com
 * @version last changed: Jun 30, 2004
 * @see
 * @since Jun 30, 2004
 */
class WarpPathWindow(path: WarpPath, searchRadius: Int) :
    SearchWindow(path.get(path.size() - 1).col + 1, path.get(path.size() - 1).row + 1) {
    companion object {
        // CONSTANTS
        private const val defaultRadius = 0
    }

    // CONSTRUCTORS
    init {
        for (p in 0 until path.size()) super.markVisited(path.get(p).col, path.get(p).row)
        super.expandWindow(searchRadius)
    } // end Constructor
} // end class WarpPathWindow


class WarpPath() {
    // DATA
    private val tsIindexes // ArrayList of Integer
            : ArrayList<Any?>
    private val tsJindexes // ArrayList of Integer
            : ArrayList<Any?>

    constructor(initialCapacity: Int) : this() {
        tsIindexes.ensureCapacity(initialCapacity)
        tsJindexes.ensureCapacity(initialCapacity)
    }

    constructor(inputFile: String) : this() {
        try {
            // Record the Label names (fropm the top row.of the input file).
            val br =
                BufferedReader(FileReader(inputFile)) // open the input file

            // Read Cluster assignments.
            var line: String?
            while (br.readLine().also { line = it } != null) // read lines until end of file
            {
                val st = StringTokenizer(line, ",", false)
                if (st.countTokens() == 2) {
                    tsIindexes.add(st.nextToken())
                    tsJindexes.add(st.nextToken())
                } else throw InternalError(
                    """
                        The Warp Path File '$inputFile' has an incorrect format.  There must be
                        two numbers per line separated by commas
                        """.trimIndent()
                )
            } // end while loop
        } catch (e: FileNotFoundException) {
            throw InternalError("ERROR:  The file '$inputFile' was not found.")
        } catch (e: IOException) {
            throw InternalError("ERROR:  Problem reading the file '$inputFile'.")
        } // end try block
    }

    // FUNCTIONS
    fun size(): Int {
        return tsIindexes.size
    }

    fun minI(): Int {
        return (tsIindexes[0] as Int).toInt()
    }

    fun minJ(): Int {
        return (tsJindexes[0] as Int).toInt()
    }

    fun maxI(): Int {
        return (tsIindexes[tsIindexes.size - 1] as Int).toInt()
    }

    fun maxJ(): Int {
        return (tsJindexes[tsJindexes.size - 1] as Int).toInt()
    }

    fun addFirst(i: Int, j: Int) {
        tsIindexes.add(0, i)
        tsJindexes.add(0, j)
    }

    fun addLast(i: Int, j: Int) {
        tsIindexes.add(i)
        tsJindexes.add(j)
    }

    fun getMatchingIndexesForI(i: Int): ArrayList<*> {
        var index = tsIindexes.indexOf(i)
        if (index < 0) throw InternalError(
            "ERROR:  index '" + i + " is not in the " +
                    "warp path."
        )
        val matchingJs: ArrayList<Any?> = ArrayList<Any?>()
        while (index < tsIindexes.size && tsIindexes[index] == i) matchingJs.add(
            tsJindexes[index++]
        )
        return matchingJs
    } // end getMatchingIndexesForI(int i)

    fun getMatchingIndexesForJ(j: Int): ArrayList<*> {
        var index = tsJindexes.indexOf(j)
        if (index < 0) throw InternalError(
            "ERROR:  index '" + j + " is not in the " +
                    "warp path."
        )
        val matchingIs: ArrayList<Any?> = ArrayList<Any?>()
        while (index < tsJindexes.size && tsJindexes[index] == j) matchingIs.add(
            tsIindexes[index++]
        )
        return matchingIs
    } // end getMatchingIndexesForI(int i)

    // Create a new WarpPath that is the same as THIS WarpPath, but J is the reference template, rather than I.
    fun invertedCopy(): WarpPath {
        val newWarpPath = WarpPath()
        for (x in tsIindexes.indices) newWarpPath.addLast(
            (tsJindexes[x] as Int).toInt(),
            (tsIindexes[x] as Int).toInt()
        )
        return newWarpPath
    }

    // Swap I and J so that the warp path now indicates that J is the template rather than I.
    fun invert() {
        for (x in tsIindexes.indices) {
            val temp = tsIindexes[x]
            tsIindexes.set(x, tsJindexes[x])
            tsJindexes.set(x, temp)
        }
    } // end invert()

    operator fun get(index: Int): ColMajorCell {
        return if (index > size() || index < 0) throw NoSuchElementException() else ColMajorCell(
            (tsIindexes[index] as Int).toInt(),
            (tsJindexes[index] as Int).toInt()
        )
    }

    override fun toString(): String {
        val outStr = StringBuffer("[")
        for (x in tsIindexes.indices) {
            outStr.append("(" + tsIindexes[x] + "," + tsJindexes[x] + ")")
            if (x < tsIindexes.size - 1) outStr.append(",")
        } // end for loop
        return String(outStr.append("]"))
    } // end toString()

    override fun equals(obj: Any?): Boolean {
        return if (obj is WarpPath) // trivial false test
        {
            val p = obj
            if (p.size() == size() && p.maxI() == maxI() && p.maxJ() == maxJ()) // less trivial reject
            {
                // Compare each value in the warp path for equality
                for (x in 0 until size()) if (!this[x].equals(p[x])) return false
                true
            } else false
        } else false
    } // end equals

    override fun hashCode(): Int {
        return tsIindexes.hashCode() * tsJindexes.hashCode()
    }

    // CONSTRUCTORS
    init {
        tsIindexes = ArrayList<Any?>()
        tsJindexes = ArrayList<Any?>()
    }
} // end class WarpPath


class PAA(ts: TimeSeries, shrunkSize: Int) : TimeSeries() {
    // PRIVATE DATA
    private val aggPtSize // ArrayList of Integer
            : IntArray
    private val originalLength: Int
    fun originalSize(): Int {
        return originalLength
    }

    fun aggregatePtSize(ptIndex: Int): Int {
        return aggPtSize[ptIndex]
    }

    override fun toString(): String {
        return """(${originalLength} point time series represented as ${this.size()} points) ${super.toString()}"""
    } // end toString()

    init {
        if (shrunkSize > ts.size()) throw InternalError(
            """ERROR:  The size of an aggregate representation may not be largerr than the 
original time series (shrunkSize=$shrunkSize , origSize=${ts.size()})."""
        )
        if (shrunkSize <= 0) throw InternalError(
            """
                ERROR:  The size of an aggregate representation must be greater than zero and 
                no larger than the original time series.
                """.trimIndent()
        )

        // Initialize private data.
        originalLength = ts.size()
        aggPtSize = IntArray(shrunkSize)

        // Ensures that the data structure storing the time series will not need
        //    to be expanded more than once.  (not necessary, for optimization)
        super.setMaxCapacity(shrunkSize)

        // Initialize the new aggregate time series.
        this.setLabels(ts.getLabels())

        // Determine the size of each sampled point. (may be a fraction)
        val reducedPtSize = ts.size() as Double / shrunkSize.toDouble()

        // Variables that keep track of the range of points being averaged into a single point.
        var ptToReadFrom = 0
        var ptToReadTo: Int


        // Keep averaging ranges of points into aggregate points until all of the data is averaged.
        while (ptToReadFrom < ts.size().toInt()) {
            ptToReadTo =
                Math.round(reducedPtSize * (this.size() + 1))
                    .toInt() - 1 // determine end of current range
            val ptsToRead = ptToReadTo - ptToReadFrom + 1

            // Keep track of the sum of all the values being averaged to create a single point.
            var timeSum = 0.0
            val measurementSums = DoubleArray(ts.numOfDimensions())

            // Sum all of the values over the range ptToReadFrom...ptToReadFrom.
            for (pt in ptToReadFrom..ptToReadTo) {
                val currentPoint: DoubleArray = ts.getMeasurementVector(pt)
                timeSum += ts.getTimeAtNthPoint(pt)
                for (dim in 0 until ts.numOfDimensions()) measurementSums[dim] += currentPoint[dim]
            } // end for loop

            // Determine the average value over the range ptToReadFrom...ptToReadFrom.
            timeSum = timeSum / ptsToRead
            for (dim in 0 until ts.numOfDimensions()) measurementSums[dim] =
                measurementSums[dim] / ptsToRead // find the average of each measurement

            // Add the computed average value to the aggregate approximation.
            aggPtSize[super.size()] = ptsToRead
            this.addLast(timeSum, TimeSeriesPoint(measurementSums))
            ptToReadFrom =
                ptToReadTo + 1 // next window of points to average startw where the last window ended
        } // end while loop
    } // end Constructor
} // end class PAA


object DTW {
    // FUNCTIONS
    fun calcWarpCost(
        path: WarpPath,
        tsI: TimeSeries,
        tsJ: TimeSeries,
        distFn: DistanceFunction
    ): Double {
        var totalCost = 0.0
        for (p in 0 until path.size()) {
            val currWarp: ColMajorCell = path.get(p)
            totalCost += distFn.calcDistance(
                tsI.getMeasurementVector(currWarp.col),
                tsJ.getMeasurementVector(currWarp.row)
            )
        }
        return totalCost
    }

    // Dynamic Time Warping where the warp path is not needed, an alternate implementation can be used that does not
    //    require the entire cost matrix to be filled and only needs 2 columns to be stored at any one time.
    fun getWarpDistBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        distFn: DistanceFunction
    ): Double {
        // The space complexity is 2*tsJ.size().  Dynamic time warping is symmetric so switching the two time series
        //    parameters does not effect the final warp cost but can reduce the space complexity by allowing tsJ to
        //    be set as the shorter time series and only requiring 2 columns of size |tsJ| rather than 2 larger columns of
        //    size |tsI|.
        if (tsI.size() < tsJ.size()) return getWarpDistBetween(tsJ, tsI, distFn)
        var lastCol = DoubleArray(tsJ.size())
        var currCol = DoubleArray(tsJ.size())
        val maxI = tsI.size() - 1
        val maxJ = tsJ.size() - 1

        // Calculate the values for the first column, from the bottom up.
        currCol[0] = distFn.calcDistance(
            tsI.getMeasurementVector(0),
            tsJ.getMeasurementVector(0)
        ) // first cell
        for (j in 1..maxJ)  // the rest of the first column
            currCol[j] = currCol[j - 1] + distFn.calcDistance(
                tsI.getMeasurementVector(0),
                tsJ.getMeasurementVector(j)
            )
        for (i in 1..maxI)  // i = columns
        {
            // Swap the references between the two arrays.
            val temp = lastCol
            lastCol = currCol
            currCol = temp

            // Calculate the value for the bottom row of the current column
            //    (i,0) = LocalCost(i,0) + GlobalCost(i-1,0)
            currCol[0] = lastCol[0] + distFn.calcDistance(
                tsI.getMeasurementVector(i),
                tsJ.getMeasurementVector(0)
            )
            for (j in 1..maxJ)  // j = rows
            {
                // (i,j) = LocalCost(i,j) + minGlobalCost{(i-1,j),(i-1,j-1),(i,j-1)}
                val minGlobalCost = Math.min(
                    lastCol[j],
                    Math.min(lastCol[j - 1], currCol[j - 1])
                )
                currCol[j] = minGlobalCost + distFn.calcDistance(
                    tsI.getMeasurementVector(i),
                    tsJ.getMeasurementVector(j)
                )
            } // end for loop
        } // end for loop

        // Minimum Cost is at (maxI,maxJ)
        return currCol[maxJ]
    } // end getWarpDistBetween(..)

    fun getWarpPathBetween(tsI: TimeSeries, tsJ: TimeSeries, distFn: DistanceFunction): WarpPath {
        return DynamicTimeWarp(tsI, tsJ, distFn).path
    }

    fun getWarpInfoBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        distFn: DistanceFunction
    ): TimeWarpInfo {
        return DynamicTimeWarp(tsI, tsJ, distFn)
    }

    private fun DynamicTimeWarp(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        distFn: DistanceFunction
    ): TimeWarpInfo {
        //     COST MATRIX:
        //   5|_|_|_|_|_|_|E| E = min Global Cost
        //   4|_|_|_|_|_|_|_| S = Start point
        //   3|_|_|_|_|_|_|_| each cell = min global cost to get to that point
        // j 2|_|_|_|_|_|_|_|
        //   1|_|_|_|_|_|_|_|
        //   0|S|_|_|_|_|_|_|
        //     0 1 2 3 4 5 6
        //            i
        //   access is M(i,j)... column-row
        val costMatrix =
            Array(tsI.size()) { DoubleArray(tsJ.size()) }
        val maxI = tsI.size() - 1
        val maxJ = tsJ.size() - 1

        // Calculate the values for the first column, from the bottom up.
        costMatrix[0][0] = distFn.calcDistance(
            tsI.getMeasurementVector(0),
            tsJ.getMeasurementVector(0)
        )
        for (j in 1..maxJ) costMatrix[0][j] =
            costMatrix[0][j - 1] + distFn.calcDistance(
                tsI.getMeasurementVector(0),
                tsJ.getMeasurementVector(j)
            )
        for (i in 1..maxI)  // i = columns
        {
            // Calculate the value for the bottom row of the current column
            //    (i,0) = LocalCost(i,0) + GlobalCost(i-1,0)
            costMatrix[i][0] = costMatrix[i - 1][0] + distFn.calcDistance(
                tsI.getMeasurementVector(i),
                tsJ.getMeasurementVector(0)
            )
            for (j in 1..maxJ)  // j = rows
            {
                // (i,j) = LocalCost(i,j) + minGlobalCost{(i-1,j),(i-1,j-1),(i,j-1)}
                val minGlobalCost = Math.min(
                    costMatrix[i - 1][j],
                    Math.min(
                        costMatrix[i - 1][j - 1],
                        costMatrix[i][j - 1]
                    )
                )
                costMatrix[i][j] = minGlobalCost + distFn.calcDistance(
                    tsI.getMeasurementVector(i),
                    tsJ.getMeasurementVector(j)
                )
            } // end for loop
        } // end for loop

/*
// writes a section of the cost matrix to a file
try
{
final PrintWriter out = new PrintWriter(new FileWriter("matrix1.csv"));
for (int j=maxJ; j>=0; j--)
{
   for (int i=0; i<=maxI; i++)
   {
      out.print(costMatrix[i][j]);
      if (i != maxI)
         out.print(",");
   }
   out.println();
}
out.flush();
out.close();
}
catch (Exception e)
{
   System.out.println(e);
   e.printStackTrace();
}
 */
        // Minimum Cost is at (maxIi,maxJ)
        val minimumCost = costMatrix[maxI][maxJ]


        // Find the Warp Path by searching the matrix from the solution at
        //    (maxI, maxJ) to the beginning at (0,0).  At each step move through
        //    the matrix 1 step left, down, or diagonal, whichever has the
        //    smallest cost.  Favoer diagonal moves and moves towards the i==j
        //    axis to break ties.
        val minCostPath = WarpPath(maxI + maxJ - 1)
        var i = maxI
        var j = maxJ
        minCostPath.addFirst(i, j)
        while (i > 0 || j > 0) {
            // Find the costs of moving in all three possible directions (left,
            //    down, and diagonal (down and left at the same time).
            val diagCost: Double
            val leftCost: Double
            val downCost: Double
            diagCost = if (i > 0 && j > 0) costMatrix[i - 1][j - 1] else Double.POSITIVE_INFINITY
            leftCost = if (i > 0) costMatrix[i - 1][j] else Double.POSITIVE_INFINITY
            downCost = if (j > 0) costMatrix[i][j - 1] else Double.POSITIVE_INFINITY

            // Determine which direction to move in.  Prefer moving diagonally and
            //    moving towards the i==j axis of the matrix if there are ties.
            if (diagCost <= leftCost && diagCost <= downCost) {
                i--
                j--
            } else if (leftCost < diagCost && leftCost < downCost) i-- else if (downCost < diagCost && downCost < leftCost) j-- else if (i <= j) // leftCost==rightCost > diagCost
                j-- else  // leftCost==rightCost > diagCost
                i--

            // Add the current step to the warp path.
            minCostPath.addFirst(i, j)
        } // end while loop
        return TimeWarpInfo(minimumCost, minCostPath)
    } // end DynamicTimeWarp(..)

    fun getWarpDistBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        window: SearchWindow,
        distFn: DistanceFunction
    ): Double {
        //     COST MATRIX:
        //   5|_|_|_|_|_|_|E| E = min Global Cost
        //   4|_|_|_|_|_|_|_| S = Start point
        //   3|_|_|_|_|_|_|_| each cell = min global cost to get to that point
        // j 2|_|_|_|_|_|_|_|
        //   1|_|_|_|_|_|_|_|
        //   0|S|_|_|_|_|_|_|
        //     0 1 2 3 4 5 6
        //            i
        //   access is M(i,j)... column-row
        val costMatrix: CostMatrix = PartialWindowMatrix(window)
        val maxI = tsI.size() - 1
        val maxJ = tsJ.size() - 1

        // Get an iterator that traverses the window cells in the order that the cost matrix is filled.
        //    (first to last row (1..maxI), bottom to top (1..MaxJ)
        val matrixIterator: Iterator<*> = window.iterator()
        while (matrixIterator.hasNext()) {
            val currentCell =
                matrixIterator.next() as ColMajorCell // current cell being filled
            val i = currentCell.col
            val j = currentCell.row
            if (i == 0 && j == 0) // bottom left cell (first row AND first column)
                costMatrix.put(
                    i,
                    j,
                    distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(0))
                ) else if (i == 0) // first column
            {
                costMatrix.put(
                    i,
                    j,
                    distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)) +
                            costMatrix[i, j - 1]
                )
            } else if (j == 0) // first row
            {
                costMatrix.put(
                    i,
                    j,
                    distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) +
                            costMatrix[i - 1, j]
                )
            } else  // not first column or first row
            {
                val minGlobalCost = Math.min(
                    costMatrix[i - 1, j],
                    Math.min(
                        costMatrix[i - 1, j - 1],
                        costMatrix[i, j - 1]
                    )
                )
                costMatrix.put(
                    i, j, minGlobalCost + distFn.calcDistance(
                        tsI.getMeasurementVector(i),
                        tsJ.getMeasurementVector(j)
                    )
                )
            } // end if
        } // end while loop

        // Minimum Cost is at (maxI, maxJ)
        return costMatrix[maxI, maxJ]
    } // end getWarpDistBetween(...)

    fun getWarpPathBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        window: SearchWindow,
        distFn: DistanceFunction
    ): WarpPath {
        return constrainedTimeWarp(tsI, tsJ, window, distFn).path
    }

    fun getWarpInfoBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        window: SearchWindow,
        distFn: DistanceFunction
    ): TimeWarpInfo {
        return constrainedTimeWarp(tsI, tsJ, window, distFn)
    }

    private fun constrainedTimeWarp(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        window: SearchWindow,
        distFn: DistanceFunction
    ): TimeWarpInfo {
        //     COST MATRIX:
        //   5|_|_|_|_|_|_|E| E = min Global Cost
        //   4|_|_|_|_|_|_|_| S = Start point
        //   3|_|_|_|_|_|_|_| each cell = min global cost to get to that point
        // j 2|_|_|_|_|_|_|_|
        //   1|_|_|_|_|_|_|_|
        //   0|S|_|_|_|_|_|_|
        //     0 1 2 3 4 5 6
        //            i
        //   access is M(i,j)... column-row
        val costMatrix = WindowMatrix(window)
        val maxI = tsI.size() - 1
        val maxJ = tsJ.size() - 1

        // Get an iterator that traverses the window cells in the order that the cost matrix is filled.
        //    (first to last row (1..maxI), bottom to top (1..MaxJ)
        val matrixIterator: Iterator<*> = window.iterator()
        while (matrixIterator.hasNext()) {
            val currentCell =
                matrixIterator.next() as ColMajorCell // current cell being filled
            val i = currentCell.col
            val j = currentCell.row
            if (i == 0 && j == 0) // bottom left cell (first row AND first column)
                costMatrix.put(
                    i,
                    j,
                    distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(0))
                ) else if (i == 0) // first column
            {
                costMatrix.put(
                    i,
                    j,
                    distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)) +
                            costMatrix.get(i, j - 1)
                )
            } else if (j == 0) // first row
            {
                costMatrix.put(
                    i,
                    j,
                    distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) +
                            costMatrix.get(i - 1, j)
                )
            } else  // not first column or first row
            {
                val minGlobalCost: Double = Math.min(
                    costMatrix.get(i - 1, j),
                    Math.min(
                        costMatrix.get(i - 1, j - 1),
                        costMatrix.get(i, j - 1)
                    )
                )
                costMatrix.put(
                    i, j, minGlobalCost + distFn.calcDistance(
                        tsI.getMeasurementVector(i),
                        tsJ.getMeasurementVector(j)
                    )
                )
            } // end if
        } // end while loop

        // Minimum Cost is at (maxI, maxJ)
        val minimumCost: Double = costMatrix.get(maxI, maxJ)

/*
try
{
final PrintWriter out = new PrintWriter(new FileWriter("matrix2.csv"));
for (int j=maxJ; j>=0; j--)
{
   for (int i=0; i<=maxI; i++)
   {
      out.print(costMatrix.get(i, j));
      if (i != maxI)
         out.print(",");
   }
   out.println();
}
out.flush();
out.close();
}
catch (Exception e)
{
   System.out.println(e);
   e.printStackTrace();
}
*/
        // Find the Warp Path by searching the matrix from the solution at
        //    (maxI, maxJ) to the beginning at (0,0).  At each step move through
        //    the matrix 1 step left, down, or diagonal, whichever has the
        //    smallest cost.  Favoer diagonal moves and moves towards the i==j
        //    axis to break ties.
        val minCostPath = WarpPath(maxI + maxJ - 1)
        var i = maxI
        var j = maxJ
        minCostPath.addFirst(i, j)
        while (i > 0 || j > 0) {
            // Find the costs of moving in all three possible directions (left,
            //    down, and diagonal (down and left at the same time).
            val diagCost: Double
            val leftCost: Double
            val downCost: Double
            diagCost =
                if (i > 0 && j > 0) costMatrix.get(i - 1, j - 1) else Double.POSITIVE_INFINITY
            leftCost = if (i > 0) costMatrix.get(i - 1, j) else Double.POSITIVE_INFINITY
            downCost = if (j > 0) costMatrix.get(i, j - 1) else Double.POSITIVE_INFINITY

            // Determine which direction to move in.  Prefer moving diagonally and
            //    moving towards the i==j axis of the matrix if there are ties.
            if (diagCost <= leftCost && diagCost <= downCost) {
                i--
                j--
            } else if (leftCost < diagCost && leftCost < downCost) i-- else if (downCost < diagCost && downCost < leftCost) j-- else if (i <= j) // leftCost==rightCost > diagCost
                j-- else  // leftCost==rightCost > diagCost
                i--

            // Add the current step to the warp path.
            minCostPath.addFirst(i, j)
        } // end while loop

        // Free any rescources associated with the costMatrix (a swap file may have been created if the swa file did not
        //    fit into main memory).
        costMatrix.freeMem()
        return TimeWarpInfo(minimumCost, minCostPath)
    } // end ConstrainedTimeWarp
} // end class DtwTest


object FastDTW {
    // CONSTANTS
    const val DEFAULT_SEARCH_RADIUS = 1
    fun getWarpDistBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        distFn: DistanceFunction
    ): Double {
        return fastDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn).distance
    }

    fun getWarpDistBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        searchRadius: Int,
        distFn: DistanceFunction
    ): Double {
        return fastDTW(tsI, tsJ, searchRadius, distFn).distance
    }

    fun getWarpPathBetween(tsI: TimeSeries, tsJ: TimeSeries, distFn: DistanceFunction): WarpPath {
        return fastDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn).path
    }

    fun getWarpPathBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        searchRadius: Int,
        distFn: DistanceFunction
    ): WarpPath {
        return fastDTW(tsI, tsJ, searchRadius, distFn).path
    }

    fun getWarpInfoBetween(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        searchRadius: Int,
        distFn: DistanceFunction
    ): TimeWarpInfo {
        return fastDTW(tsI, tsJ, searchRadius, distFn)
    }

    private fun fastDTW(
        tsI: TimeSeries,
        tsJ: TimeSeries,
        searchRadius: Int,
        distFn: DistanceFunction
    ): TimeWarpInfo {
        var searchRadius = searchRadius
        if (searchRadius < 0) searchRadius = 0
        val minTSsize = searchRadius + 2
        return if (tsI.size() <= minTSsize || tsJ.size() <= minTSsize) {
            // Perform full Dynamic Time Warping.
            getWarpInfoBetween(tsI, tsJ, distFn)
        } else {
            val resolutionFactor = 2.0
            val shrunkI = PAA(tsI, (tsI.size() / resolutionFactor).toInt())
            val shrunkJ = PAA(tsJ, (tsJ.size() / resolutionFactor).toInt())

            // Determine the search window that constrains the area of the cost matrix that will be evaluated based on
            //    the warp path found at the previous resolution (smaller time series).
            val window: SearchWindow = ExpandedResWindow(
                tsI, tsJ, shrunkI, shrunkJ,
                getWarpPathBetween(shrunkI, shrunkJ, searchRadius, distFn),
                searchRadius
            )

            // Find the optimal warp path through this search window constraint.
            getWarpInfoBetween(tsI, tsJ, window, distFn)
        } // end if
    } // end recFastDTW(...)
} // end class fastDTW
