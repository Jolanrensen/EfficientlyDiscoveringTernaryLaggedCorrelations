package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations

import kotlin.properties.ReadWriteProperty
import kotlin.reflect.KProperty

/**
 * Can be used to create a delegatable reference to, for instance, a primitive type.
 * */
class Delegate<T>(var value: T) : ReadWriteProperty<Any?, T> {
    override fun getValue(thisRef: Any?, property: KProperty<*>): T = value

    override fun setValue(thisRef: Any?, property: KProperty<*>, value: T) {
        this.value = value
    }
}

inline fun <reified T> delegateOf(initial: T): Delegate<T> = Delegate(initial)

