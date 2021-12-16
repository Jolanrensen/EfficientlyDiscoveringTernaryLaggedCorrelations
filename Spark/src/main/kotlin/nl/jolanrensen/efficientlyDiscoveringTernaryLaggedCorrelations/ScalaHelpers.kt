package nl.jolanrensen.thesis

import scala.*
import scala.collection.TraversableOnce
import scala.collection.generic.Growable

operator fun <A> Growable<A>.plusAssign(elem: A) {
    `$plus$eq`(elem)
}

operator fun <A> Growable<A>.plusAssign(xs: TraversableOnce<A>) {
    `$plus$plus$eq`(xs)
}

inline fun <S, T> Tuple1<S>.map(block: (S) -> T): Tuple1<T> = Tuple1<T>(block(_1()))
inline fun <S, T> Tuple2<S, S>.map(block: (S) -> T): Tuple2<T, T> = Tuple2<T, T>(block(_1()), block(_2()))
inline fun <S, T> Tuple3<S, S, S>.map(block: (S) -> T): Tuple3<T, T, T> = Tuple3<T, T, T>(block(_1()), block(_2()), block(_3()))
inline fun <S, T> Tuple4<S, S, S, S>.map(block: (S) -> T): Tuple4<T, T, T, T> = Tuple4<T, T, T, T>(block(_1()), block(_2()), block(_3()), block(_4()))
inline fun <S, T> Tuple5<S, S, S, S, S>.map(block: (S) -> T): Tuple5<T, T, T, T, T> = Tuple5<T, T, T, T, T>(block(_1()), block(_2()), block(_3()), block(_4()), block(_5()))
// TODO maybe add others